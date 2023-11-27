#!/usr/bin/perl

use forks;
use forks::shared;

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil :sys_wait_h);
use Compress::Zlib;
use File::Copy;
use File::Spec;
use File::Which;
use File::Basename;
use File::Temp qw(tempfile);
use IPC::Shareable qw(:lock);
use Storable qw(nfreeze thaw);

BEGIN {
    binmode(STDIN);
    binmode(STDOUT);
    select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
    select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

    $SIG{INT} = sub {exit(130)};

    mkdir("$ENV{HOME}/.Inline") if(! -d "$ENV{HOME}/.Inline");
}

#use Inline (Config => DIRECTORY => File::Spec->tmpdir());
use Inline qw(C);

my ($exe) = $0 =~ /([^\/]+)$/;
my $usage = "
Usage:

     $exe <bam_file>

     Organizes reads from existing BAM file into pairs.

Options:
     cpus|c     <INT>     CPUs to use for conversion

     level      <INT>     Compression level for output

     out       <PATH>     Output file (STDOUT otherwise)

     sam                  Output in SAM format

     help|?               Prints this usage statement

";

#get options from command line
my $outfile;
my $sam;
my $level;
my $cpus = 1;
GetOptions("cpus|c=i" => \$cpus,
	   "out=s"    => \$outfile,
	   "level=i"  => \$level,
	   "sam"      => \$sam,
           "help|?"   => sub{print $usage; exit(0)});
$level = -1 if($sam && !$level); #no blocks on SAM
$level = 4 if(! defined($level));

my $file = $ARGV[0];
if(!$file){
    print $usage;
    exit(0);
}
die "ERROR: The file $file does not exist\n" if(! -f $file);

#make sections of right size
my $size = (stat($file))[7];
my $count = ceil($size/33554432); #32Mb sections
my $window = ceil($size/$count);
my @temp_files : shared; #holds intermediate data from threads 
my @sections = (1..$count);
share(@sections); #share is non destructive for forks::shared

#prepare output file(s)
my $osize = 1.5*$size;
if($outfile){ #preallocate
    get_handle('>', $outfile); #initialize
    close_handles(); #clear handle
    truncate($outfile, $osize);
}

#make shared values I need
our $IPC_SHARE;
our $IPCHANDLE = tie($IPC_SHARE, 'IPC::Shareable', undef, { destroy => 1 });
$IPC_SHARE = 0; #output buffer offset

#get header
my $header = get_header($file);
my $header_data = ($sam) ? $header->{text} : $header->{data};
print_async(zipit($level,\$header_data), $outfile, 0);

#run all threads on data
my %param = (SIZE => $size, COUNT => $count, WINDOW => $window, SECTIONS => \@sections,
	     INFILE => $file, OUTFILE => $outfile, LEVEL => $level, HEADER => $header);
for(my $i = 1; $i < $cpus; $i++){
    $param{ID} = $i;
    threads->new({'context' => 'scalar'}, \&thread_run, \%param);
}
$param{ID} = 0;
thread_run(\%param); #let main process join in

#clean up threads
$_->join foreach(threads->list());

#finish output
print_async(\ (eof_block()), $outfile, 0);
close_handles(); #clear handles
truncate($outfile, $IPC_SHARE) if($outfile && -f $outfile); #make file proper size

#finished
exit(0);

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------

#get header from bam
sub get_header {
    my $file = shift;
    my $validate = 1; #just leave it on

    #declare variables
    my $header;
    my $header_off;
    my $header_block;
    
    #open input file
    my $IN = get_handle('<', $file);
    binmode($IN);
    
    #get header
    inflate(readblock($IN), \$header_block, $validate);
    
    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
	if($magic ne "BAM\1");
    
    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
	inflate(readblock($IN), \$header_block, $validate);
    }
    my $text = substr($header_block, $header_off, $l_text);
    $header_off += $l_text;
    $text =~ s/\0$//g; #remove null padding
    $header->{text} = $text;
    
    my $n_ref = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    $header->{n_ref} = $n_ref;

    #get each reference
    for(my $i = 0; $i < $n_ref; $i++){
	while(length($header_block)-$header_off < 4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $l_name = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	
	while(length($header_block)-$header_off < $l_name+4){ #grow if needed
	    inflate(readblock($IN), \$header_block, $validate);
	}
	my $name = substr($header_block, $header_off, $l_name);
	$header_off += $l_name;
	$name =~ s/\0$//g; #remove null padding
	$header->{ref}[$i]{name} = $name;
	$header->{refid2name}[$i] = $name;

	my $l_ref = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	$header->{ref}[$i]{l_ref} = $l_ref;
    }
    $header->{data} = substr($header_block,0,$header_off,''); #chop off header from block
    $header_off = tell($IN); #make offset be position of first block after header

    $header->{offset} = $header_off; #offset of first block imediately following the header
    $header->{block}  = $header_block; #piece of alignment accidentally stuck inside header block

    return $header;
}

#convenience method to see if I am in a thread or not
{my $tid; BEGIN{$tid = threads->tid if exists $INC{'threads.pm'}};
sub is_thread {
    my $is_thread = 1 if(defined($tid) && $tid != threads->tid);
    return $is_thread;
}}

#main fastq conversion done here by each process
sub thread_run {
    my $param = shift;

    #fix thread signalling
    if(is_thread){
	binmode(STDIN);
	binmode(STDOUT);
	select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
	select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

	$SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
	$SIG{INT} = sub {exit(130)};
    }

    #load parameters
    my $sections = $param->{SECTIONS};
    my $size     = $param->{SIZE};
    my $count    = $param->{COUNT};
    my $window   = $param->{WINDOW};
    my $file     = $param->{INFILE};
    my $outfile  = $param->{OUTFILE};
    my $level    = $param->{LEVEL};
    my $header   = $param->{HEADER};
    my $id       = $param->{ID};
    my $validate = 1; #just leave this on

    #declare variables
    my @buffer = map { {} } (1..$cpus); #paired read buffer (array of hash refs)
    my $bam_ref;
    {my $string = ''; $bam_ref = \$string;} #set reference
    my $bgz_ref;
    {my $string = ''; $bgz_ref = \$string;} #set reference
    
    #open input file
    close_handles(); #clear existing handles to avoid duplication among threads
    my $IN = get_handle('<', $file);
    binmode($IN);

    #now read alignment sections
    while(my $section = shift @$sections){
	my $start = ($section-1) * $window;
	my $end   = $section * $window;
	$start = $header->{offset} if($start < $header->{offset});
	$end = $size if($end > $size); #don't go after end

	#adjust to start of next block
	seek_next_block($IN, $start, 0);
	
	#get first block in section
	my $block = ''; #piece of alignment stuck at end of header
	$block = $header->{block} if($section == 1);
	inflate(readblock($IN), \$block, $validate);

	#find first alignment
	my $offset = tell_next_alignment(\$block, $header);
	while(! defined($offset) && tell($IN) < $start+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$block, $validate);
	    $offset = tell_next_alignment(\$block, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($block));
	substr($block, 0, $offset) = ''; #chop off leading data

	#get all blocks in section
	while(tell($IN) < $end){
	    inflate(readblock($IN), \$block, $validate);
	    pairbam(\$block, $bam_ref, \@buffer);
	}

	#get last block and adjust split terminating line
	my $tail = '';
	inflate(readblock($IN), \$tail, $validate);
	$offset = tell_next_alignment(\$tail, $header);
	while(! defined($offset) && tell($IN) < $end+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$tail, $validate);
	    $offset = tell_next_alignment(\$tail, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($tail));
	$block .= substr($tail, 0, $offset) if($offset); #chop off trailing data

	#convert to fastq
	pairbam(\$block, $bam_ref, \@buffer);

	print_async(zipit($level, $bam_ref, $bgz_ref), $outfile, 1) if(length($$bam_ref) >= 33554432);
    }

    #thread stores buffer to avoid high memory usage
    my ($tfh, $tfile) = tempfile("pairbam\_$id\_XXXXX", CLEANUP => 1, DIR => File::Spec->tmpdir);
    for(my $i = 0; $i < @buffer; $i++){
	next if(! keys %{$buffer[$i]}); #skip empty hashes
	$buffer[$i] = nfreeze($buffer[$i]);
	print $tfh pack('S<L<', $i, length($buffer[$i])); #(fasthashid, length)
	print $tfh $buffer[$i];
	undef $buffer[$i];
    }
    close($tfh);

    #share file among threads
    push(@temp_files, $tfile);
    sleep 1 while(@temp_files < $cpus);

    #gather and try and to colapse remaining temp_file buffers
    foreach my $file (@temp_files){
	open($tfh, "< $file");
	my $meta;
	while(read($tfh, $meta, 6)){
	    my ($fhid, $len) = unpack('S<L<', $meta); #is this the buffer section I want
	    if($fhid % $cpus != $id){ #skip past unwanted buffer sections
		seek($tfh, $len, 1);
		next;
	    }

	    #retrieve buffer section
	    my $hash2;
	    read($tfh, $hash2, $len);
	    $hash2 = thaw($hash2);
	    my $hash1 = $buffer[$fhid];
	    if(!$hash1){
		$buffer[$fhid] = $hash2;
		next;
	    }
	    
	    #merge in retrieved buffer section
	    foreach my $read_name (keys %$hash2){
		#initialize read data if not in hash yet
		my $read_data2 = $hash2->{$read_name};
		if(!exists($hash1->{$read_name})){ #copy over to local hash
		    $hash1->{$read_name} = $read_data2;
		    next;
		}
		
		#merge read data into existing hash
		my $complete = 1; #flag if all data is here
		my $read_data1 = $hash1->{$read_name}; #(pair_count, pair1, pair2, ...)
		for(my $i = 1; $i <= $read_data1->[0]; $i++){
		    #merge read_data2 into read_data1
		    if($read_data2->[$i]){
			if(!$read_data1->[$i]){ #just copy it over
			    $read_data1->[$i] = $read_data2->[$i];
			}
			else{ #merge pieces
			    if($read_data2->[$i][0] == 1){ #only override if reset to 1
				$read_data1->[$i][0] = $read_data2->[$i][0]; #count of expected alignments for pair
			    }
			    $read_data1->[$i][1] += $read_data2->[$i][1]; #count found thus far for pair
			    $read_data1->[$i][2] .= $read_data2->[$i][2]; #primary alignment
			    $read_data1->[$i][3] .= $read_data2->[$i][3]; #secondary and complementary alignments
			}
		    }
		    
		    if(!$read_data1->[$i]){
			$complete = 0;
		    }
		    else{
			$complete = 0 if($read_data1->[$i][1] < $read_data1->[$i][0]); #found less than expected
		    }
		}
		
		#all data collected from all pairs
		if($complete){
		    for(my $i = 1; $i <= $read_data1->[0]; $i++){
			$$bam_ref .= $read_data1->[$i][2];
			$$bam_ref .= $read_data1->[$i][3] if($read_data1->[1] > 1);
		    }
		    delete($hash1->{$read_name}); #drop completed read
		}
	    }
	    
	    print_async(zipit($level, $bam_ref, $bgz_ref), $outfile, 1) if(length($$bam_ref) >= 33554432);
	}
	close($tfh);
    }

    #print leftovers
    print_async(zipit($level, $bam_ref, $bgz_ref), $outfile, $level, 0);

    #anything that is still in hash is either missing mates or was realigned
    for (my $i = 0; $i < @buffer; $i++){
	next if(!$buffer[$i]);
	my $hash1 = $buffer[$i];
	foreach my $read_name (keys %$hash1){
	    warn "WARNING: Read is either missing its pair or lost information when realigned: $read_name\n";
	}
    }

    #clear filehandles
    close_handles();

    return;
}

our %HANDLES;
sub get_handle {
    my $mode = shift;
    my $file = shift || 'PIPE';

    if(!$file){
	binmode(STDIN);
	binmode(STDOUT);
	return \*STDOUT if($mode eq '>' || $mode eq '>>' || $mode eq '+<');
	return \*STDIN if($mode eq '<');
    }
    if(!$HANDLES{$mode}{$file}){
	open(my $FH, $mode, $file);
	$HANDLES{$mode}{$file} = $FH;
	binmode($FH);
    }

    return $HANDLES{$mode}{$file};
}
sub close_handles {
    close($_) foreach(map {values %$_} values %HANDLES);
    undef %HANDLES;
}

sub zipit {
    my $level = shift;
    my $iref  = shift;
    my $oref  = shift;

    #compress data before printing
    if($level >= 0){
	if(!$oref){
	    my $string = '';
	    $oref = \$string;
	}
	deflate(\ (substr($$iref, 0, 65280, '')), $oref, $level) while(length($$iref));
	return $oref;
    }

    return $iref;
}

sub print_async {
    my $ref      = shift;
    my $outfile  = shift;
    my $async    = shift;

    #lock and change file size if necessary
    if($async){
	return unless($IPCHANDLE->shlock(LOCK_EX|LOCK_NB));
    }
    else{
	sleep 0.1 while(!$IPCHANDLE->shlock(LOCK_EX));
    }

    #output first pair or interleaved pairs
    my $loc = $IPC_SHARE;
    $IPC_SHARE += length($$ref);
    if(!$outfile){
	#keep trying until it prints everything
	my $offset = 0;
	$offset += syswrite(STDOUT, $$ref, length($$ref), $offset) while(length($$ref) > $offset);
	$$ref = '';
	$IPCHANDLE->shunlock;
	
	return;
    }
    elsif(-f $outfile && $IPC_SHARE > $osize){ #grow
	$osize = $IPC_SHARE+33554432; #with extra 32MB buffer
	truncate($outfile, $osize);
    }
    $IPCHANDLE->shunlock;

    #write data once lock is no longer needed
    #output first pair or interleaved pairs
    return if(! length($$ref));
    
    my $OUT = get_handle('+<', $outfile);
    binmode($OUT);
    sysseek($OUT, $loc, 0);
    my $offset = 0;
    $offset += syswrite($OUT, $$ref, length($$ref), $offset) while(length($$ref) > $offset);
    $$ref = '';

    return;
}
    
{#predeclare values for efficienty
my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen,
    $bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,$seq,@qual,$qual);
my ($tag, $val_type, $sa_value, $val_len, $sub_type, $sub_len);
my ($cl_ref,$cl_query,$match);
my ($substr_off, $data_len, $block_offset, $align_offset, $align_length);
my ($hash_key, $read_data, $pair, $sa_count, $sa_offset, $complete);
my ($data, $bam_ref, $buffer, $buf_size);

#destructive to data in block
sub pairbam {
    $data = shift; #ref
    $bam_ref = shift; #ref
    $buffer  = shift; #ref (buffer must be preallocated)
    $buf_size = scalar(@$buffer);

    die "ERROR: Failure to preallocate buffer before calling pairbam\n" if(!$buf_size);

    #collect and process bam alignments
    $data_len = length($$data);
    $block_offset = 0;
    while($block_offset < $data_len - 4){ #continue until near end of block
	$substr_off = $block_offset;
	
	($block_size) = unpack('l<', substr($$data, $substr_off, 4));
	$substr_off += 4;
	last if($substr_off+$block_size > $data_len); #current alignment is spit across next block
	$align_offset = $block_offset;
	$align_length = $block_size+4;
	$block_offset = $substr_off+$block_size; #end of current alignmnet

	#skip past refID and pos
	$substr_off += 8;

	#get bin_mq_nl, flag_nc, and l_seq 
	($bin_mq_nl,$flag_nc,$l_seq) = unpack('L<L<l<', substr($$data, $substr_off, 12));
	$substr_off += 12;

	#skip past next_refID, next_pos, and tlen
	$substr_off += 12;

	#bin_mq_nl processing into sub values
	#$bin = ($bin_mq_nl>>16);
	#$mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1 
	$l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	
	#flag_nc processing into sub values
	$flag = ($flag_nc>>16); #flag is probably not usful here
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name);
	$substr_off += $l_read_name;
	chop($read_name); #removes trailing null

	#skip past unneeded values
	$substr_off += $n_cigar_op*4; #skip past cigar string
        $substr_off += int(($l_seq+1)/2); #skip past seq string
        $substr_off += $l_seq; #skip past quality string

	#process arbitrary attributes
	$sa_value = '';
	while($substr_off < $block_offset){
	    $tag =  substr($$data, $substr_off, 2);
	    $substr_off += 2;
	    $val_type = substr($$data, $substr_off, 1);
	    $substr_off += 1;

	    $val_len = 0;
	    if   ($val_type eq 'A'){ $val_len = 1 }
	    elsif($val_type eq 'c'){ $val_len = 1 }
	    elsif($val_type eq 'C'){ $val_len = 1 }
	    elsif($val_type eq 's'){ $val_len = 2 }
	    elsif($val_type eq 'S'){ $val_len = 2 }
	    elsif($val_type eq 'i'){ $val_len = 4 }
	    elsif($val_type eq 'I'){ $val_len = 4 }
	    elsif($val_type eq 'f'){ $val_len = 4 }
	    elsif($val_type eq 'd'){ $val_len = 8 }
	    elsif($val_type eq 'Z' || $val_type eq 'H'){
		$val_len = (index($$data,"\0",$substr_off) - $substr_off)+1; #plus 1 for null                                                                                                                                                                                                                           
	    }
	    elsif($val_type eq 'B'){
		$sub_type = substr($$data, $substr_off, 1);
		$substr_off += 1;
		$sub_len = unpack('l<', substr($$data, $substr_off, 4));
		$substr_off += 4;

		if   ($sub_type eq 'c'){ $val_len = 1*$sub_len }
		elsif($sub_type eq 'C'){ $val_len = 1*$sub_len }
		elsif($sub_type eq 's'){ $val_len = 2*$sub_len }
		elsif($sub_type eq 'S'){ $val_len = 2*$sub_len }
		elsif($sub_type eq 'i'){ $val_len = 4*$sub_len }
		elsif($sub_type eq 'I'){ $val_len = 4*$sub_len }
		elsif($sub_type eq 'f'){ $val_len = 4*$sub_len }
	    }

	    #process tag
	    if($tag eq 'SA'){
		die "ERROR: Wrong datatype for SA value\n" if($val_type ne 'Z');
		$sa_value = substr($$data, $substr_off, $val_len-1); #-1 to ignore null (no need to call chomp)
		$substr_off += $val_len;
		last; #short circuit loop
	    }
	    else{ #skip past value since it is not the right one
		$substr_off += $val_len;
	    }
	}

	#group pairs together
	$hash_key = fasthash($read_name) % $buf_size;
	$buffer->[$hash_key]{$read_name} ||= [(($flag & 1) ?  2 : 1)]; #assume 2 or 1 #temp
	$read_data = $buffer->[$hash_key]{$read_name}; #(pair_count, pair1, pair2, ...)
	
	#count of expected alignments for pair
	$pair = ($flag & 128) ? 2 : 1; #which pair am I?
	if(!$sa_value || !defined($read_data->[$pair][0])){
	    $sa_count = 0;
	    $sa_offset = 0;
	    while($sa_value){
		$sa_offset = index($sa_value, ';', $sa_offset);
		last if($sa_offset == -1);
		$sa_count++;
		$sa_offset++;
	    }
	    $read_data->[$pair] = [1+$sa_count, 0, '', '']; #(expected, found, primary, secondary)
	}
	$read_data->[$pair][1]++; #count of alignments found thus far
	
	if($flag & 2048 || $flag & 256){ #secondary alignment
	    $read_data->[$pair][3] .= substr($$data, $align_offset, $align_length);
	}
	else{ #primary alignment
	    $read_data->[$pair][2] .= substr($$data, $align_offset, $align_length);
	}
	
	#check if read is complete
	$complete = 1; #flag if all data is here
	for(my $i = 1; $i <= $read_data->[0]; $i++){
	    if(!$read_data->[$i]){ #missing pair
		$complete = 0;
                last;
	    }
	    elsif($read_data->[$i][1] < $read_data->[$i][0]){ #found less than expected
		$complete = 0;
		last;
	    }
	}
	
	#all data collected from all pairs
	if($complete){
	    for(my $i = 1; $i <= $read_data->[0]; $i++){
		$$bam_ref .= $read_data->[$i][2];
		$$bam_ref .= $read_data->[$i][3] if($read_data->[1] > 1);
	    }
	    delete($buffer->[$hash_key]{$read_name}); #drop completed read
	}
    }
    
    $$data = substr($$data, $block_offset); #remove what I processed off of the block

    return;
}}

{#predeclare values for efficienty
my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen,
    $bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,$cigar,$seq,$qual);
my ($qname, $rnext, $count);
my ($tag,$val_type,$value,$val_len,$sub_type,$sub_len);
my ($cl_ref,$cl_query,$match);
my ($substr_off,$data,$data_len,$block_offset);

#destructive to data in bam
sub bam2sam {
    my $bam_ref = shift;
    my $sam_ref = shift;
    my $refid2name = shift;

    #make new sam ref if not supplied
    if(!$sam_ref){
	my $string = '';
	$sam_ref = \$string;
    }

    #convert each alignment to sam format
    $data_len = length($$bam_ref);
    $block_offset = 0;
    while($block_offset < $data_len){ #continue until end of block
	$substr_off = $block_offset;

	($block_size) = unpack('l<', substr($$data, $substr_off, 4));
	$substr_off += 4;
	$block_offset = $substr_off+$block_size; #end of current alignmnet 

	($refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen)
	    = unpack('l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 32));
	$qname = ($refID == -1) ? '*' : $refid2name->[$refID];
	$substr_off += 32;

        #bin_mq_nl processing into sub values
        $bin = ($bin_mq_nl>>16);
        $mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1
        $l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1

	#flag_nc processing into sub values
        $flag = ($flag_nc>>16); #flag is probably not usful here
        $n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
        next if($flag & 2816); #skip secondary, supplemental, and vendor failed alignments

        #get read name
        $read_name = substr($$data, $substr_off, $l_read_name-1); #remove null
        $substr_off += $l_read_name;

	#process cigar into sub values
	$cigar = substr($$data, $substr_off, $n_cigar_op*4);
	$substr_off += $n_cigar_op*4;
	convert_cigar(\$cigar, $n_cigar_op);

	#process seq string
        $count = int(($l_seq+1)/2);
        $seq = substr($$data, $substr_off, $count);
        $substr_off += $count;
        convert_seq(\$seq, $l_seq, ($flag & 16));

        #process quality string
        $qual = substr($$data, $substr_off, $l_seq);
        $substr_off += $l_seq;
        convert_qual(\$qual, $l_seq, ($flag & 16), 0);

	#process tags
	while($substr_off < $block_offset){
	    $tag =  substr($$data, $substr_off, 2);
	    $substr_off += 2;
	    $val_type = substr($$data, $substr_off, 1);
	    $substr_off += 1;
	    
	    $val_len = 0;
	    if($val_type eq 'A'){
		$val_len = 1;
		$value = substr($$data, $substr_off, $val_len);
	    }
	    elsif($val_type eq 'Z' || $val_type eq 'H'){
		$val_len = (index($$data,"\0",$substr_off) - $substr_off)+1; #plus 1 for null
		$value = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
	    }
	    elsif($val_type eq 'B'){
		$sub_type = substr($$data, $substr_off, 1);
		$substr_off += 1;
		$sub_len = unpack('l<', substr($$data, $substr_off, 4));
		$substr_off += 4;
		
		if   ($sub_type eq 'c'){ $val_len = 1*$sub_len }
		elsif($sub_type eq 'C'){ $val_len = 1*$sub_len }
		elsif($sub_type eq 's'){ $val_len = 2*$sub_len }
		elsif($sub_type eq 'S'){ $val_len = 2*$sub_len }
		elsif($sub_type eq 'i'){ $val_len = 4*$sub_len }
		elsif($sub_type eq 'I'){ $val_len = 4*$sub_len }
		elsif($sub_type eq 'f'){ $val_len = 4*$sub_len }

		$value = join(',', unpack($sub_type.$sub_len, substr($$data, $substr_off, $val_len)));
	    }
	    else{
		if($val_type eq 'c'){ $val_len = 1 }
		elsif($val_type eq 'C'){ $val_len = 1 }
		elsif($val_type eq 's'){ $val_len = 2 }
		elsif($val_type eq 'S'){ $val_len = 2 }
		elsif($val_type eq 'i'){ $val_len = 4 }
		elsif($val_type eq 'I'){ $val_len = 4 }
		elsif($val_type eq 'f'){ $val_len = 4 }
		elsif($val_type eq 'd'){ $val_len = 8 }

		$value = unpack($val_type, substr($$data, $substr_off, $val_len));
	    }

	    #jump to next tag
	    $substr_off += $val_len;
	}
    }

    $$data = substr($$data, $block_offset); #remove what I processed off of the block

    return $sam_ref;
}}

sub seek_next_block {
    my $FH = shift; #filehandle
    my $start = shift; #where to start looking from
    my $wence = shift || 0;

    #go to position
    seek($FH, $start, $wence);
    my $pos = tell($FH);
    my $size = (stat($FH))[7]; #file size

    #not enough space left for a block
    return seek($FH, $size, 0) if($size-$start < 28);

    #read in some data to check
    my $data;
    read($FH, $data, 65536, 0); #max block

    #find block keys in data
    my $offset = 0;
    my $key1 = pack('H8', '1f8b0804'); #4 byte static header key
    my $key2 = pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes)
    while(1){
	#match first static key
	$offset = index($data, $key1, $offset);
	return seek($FH, $size, 0) if($offset == -1);
	
	#match second static key
	my $offset2 = index($data, $key2, $offset);
	return seek($FH, $size, 0) if($offset2 == -1);

	#second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
	if($offset2-$offset == 10){
	    last;
	}
	else{
	    $offset += 4; #jump over key1 match and try again
	}
    }
    $pos += $offset; #adjust file position with offset

    return seek($FH, $pos, 0);
}

#assumes file position is set to start of block
sub readblock {
    my $FH = shift;

    #read BGZF header
    my $data;
    my $stat = read($FH, $data, 18, 0);
    return undef if($stat == 0);
    die "ERROR: Failure to read BAM stream\n" if($stat != 18);

    #validate header
    my ($id1,$id2,$cm,$flg,$mtime,$xfl,$os,
	$xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $data);
    die "ERROR: Does not appear to be a BGZF file\n"
	if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
	   $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);
    
    #read compression block and footer
    my $c_data_size = $bsize-$xlen-19; #compression block size
    $stat = read($FH, $data, $c_data_size + 8, length($data)); #the +8 is for the footer 
    die "ERROR: Could not read compression block\n"
	if($stat != $c_data_size + 8);

    return \$data;
}

sub inflate {
    my $buffer = shift; #reference
    my $ref    = shift; #reference
    my $validate = shift;

    return if(!$buffer);

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 131072);
    die "ERROR: Failed to create zlib inflation object with status: $stat\n" if($stat);

    substr($$buffer, 0, 18, '');
    $$ref .= scalar($i_obj->inflate($buffer));

    #validate the length and crc32
    if($validate){
	die "ERROR: Trailing garbage in compression block\n" if(length($$buffer) != 8);

        my ($crc32, $isize) = unpack('L<L<', $$buffer);
	
        if($isize != $i_obj->total_out()){ #size does not match
            die "ERROR: The expected ISIZE of the uncompressed block does not match\n";
        }
        if($crc32 != crc32(substr($$ref, -$isize, $isize))){ #crc32 does not match
            die "ERROR: The expected CRC32 of the uncompressed block does not match\n";
        }
    }

    return;
}

sub deflate {
    my $buffer = shift; #reference
    my $ref = shift; #reference
    my $level = shift;

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 131072, -Level => $level);
    die "ERROR: Failed to create zlib deflation object with status: $stat\n" if($stat);

    my $offset = length($$ref);
    $$ref .= pack('H36', '1f8b08040000000000ff0600424302000000'); #header (BSIZE=0)
    $$ref .= scalar($d_obj->deflate($buffer)).scalar($d_obj->flush); #compressed data
    $$ref .= pack('V', crc32($buffer)); #CRC32
    $$ref .= pack('V', $d_obj->total_in()); #ISIZE
    substr($$ref, $offset+16, 2, pack('v', $d_obj->total_out+6+19)); #set final BSIZE

    return;
}

#tells the virtual offset of the next alignment record
sub tell_next_alignment {
    my $data = shift; #ref
    my $header = shift;

    return if(!$$data);

    #predeclare values for efficienty
    my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,
	$tlen,$bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,@qual);
    my ($substr_off, $l_ref, $match, $cl_ref, $cl_query);
    my $n_ref = $header->{n_ref} if($header);

    #find correct virtual offset
    for(my $i = 0; $i < length($$data)-44; $i++){
	$substr_off = $i; #shift 1 byte to the right each time through loop
	
	($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,
	 $next_refID,$next_pos,$tlen) = unpack('l<l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 36));
	$substr_off += 36;
	
	next unless(44 <= $block_size && $block_size < 65536); #reasonable assumption?
	next unless( 1 <= $l_seq && $l_seq < 10000); #reasonable assumption?
	
	#bin_mq_nl processing into sub values
	$bin = ($bin_mq_nl>>16);
	$mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1 
	$l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	next unless($bin <= 37450); #reasonable assumption?
	
	#flag_nc processing into sub values
	$flag = ($flag_nc>>16);
	next unless($flag < 4096); #only first 12 bits are used
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name);
	$substr_off += $l_read_name;
	next unless($read_name =~ /^(?:\*|[!-()+-<>-~][!-~]*)\0$/); #restricted character space
	
	#process cigar into sub values
	@cigar = unpack("L<"x$n_cigar_op, substr($$data, $substr_off, $n_cigar_op*4));
	$substr_off += $n_cigar_op*4;
	@cigar = map {[($_ & 15), ($_>>4)]} @cigar; #mask is (1<<16)-1
	$cl_ref = 0;
	$cl_query = 0;
	foreach (@cigar){
	    $match = ((0x3C1A7 >> ($_->[0] << 1)) & 3);
	    $cl_query += $_->[1] if($match & 1);
	    $cl_ref += $_->[1] if($match & 2);
	}
	next unless(!@cigar || $cl_query == $l_seq); #reasonable assumption?
	next unless(reg2bin($pos, $pos+$cl_ref) == $bin || ($pos == -1 && $bin == 0));

	#removed because TLEN is niether consistent nor calculable without both reads
	#next unless($refID != $next_refID || $pos < 0 ||
	#	    $next_pos < 0 || $next_pos - ($pos+$cl_ref) == $tlen);
	
	#additionally validate contig and positions info (necessary?)
	if($header && $header->{text}){
	    next unless(-1<= $refID && $refID < $n_ref);
	    $l_ref = ($refID != -1) ? $header->{ref}[$refID]{l_ref} : 0;
	    next unless(-1 <= $pos && $pos < $l_ref);
	    next unless(-1 <= $next_refID && $next_refID < $n_ref);
	    $l_ref = ($next_refID != -1) ? $header->{ref}[$next_refID]{l_ref} : 0;
	    next unless(-1 <= $next_pos && $next_pos < $l_ref);
	}
	
	#process seq string
	
	#process quality string
	
	#ignore aux values for now
	
	return $i; #seems to be a valid bam line
    }

    return;
}

sub eof_block {
    return pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000');
}

sub thread_validate {
    my $param = shift;

    #fix thread signalling
    if(is_thread){
        binmode(STDIN);
        binmode(STDOUT);
        select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
        select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

        $SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
        $SIG{INT} = sub {exit(130)};
    }

    #load parameters
    my $size = $param->{SIZE};
    my $count = $param->{COUNT};
    my $window = $param->{WINDOW};
    my $file = $param->{INFILE};
    my $sections = $param->{SECTIONS};
    my $validate = 1;

    #open input file
    my $IN = get_handle('<', $file);
    binmode($IN);

    #now read alignment sections
    while(my $section = shift @$sections){
        my $start = ($section-1) * $window;
        my $end   = $section * $window;
        $end = $size if($end > $size); #don't go after end

        #adjust to start of next block
        seek_next_block($IN, $start, 0) if($start > 0);

        #validate all blocks in section
        while(tell($IN) < $end){
            my $block;
            inflate(readblock($IN), \$block, $validate);
        }

        #validate one more (avoids weirdness with bad blocks at boundaries)
        if(tell($IN) < $size){
            my $block;
            inflate(readblock($IN), \$block, $validate);
        }
    }

    #clear any open file handles
    close_handles();

    return;
}

#C code starts here
 __END__
__C__

#include <stdint.h>
#include <math.h>

#define mix_fasthash(h) ({ (h) ^= (h) >> 23; (h) *= 0x2127599bf4325c37ULL; (h) ^= (h) >> 47; })

unsigned long fasthash (SV* sv) {
	char* string;
	STRLEN len;
	string = SvPV(sv, len);

	uint64_t hash = 0;
	int i;

	for(i = 0; i < len; i++){
	    hash += (uint64_t)string[i];
	    mix_fasthash(hash);
	}

	return hash;
}

#define bam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

int convert_seq(SV* ref, int len, int reverse) {
    SV* scalar = SvRV(ref);
    uint8_t* seq = (uint8_t*)SvPV_nolen(scalar);

    char* buf = (char*)malloc((len+1)*sizeof(char));
    buf[len] = 0;

    int i;
    if (reverse != 0) {
	for (i = 0; i < len; ++i){
	    buf[len - i - 1] = bam_nt16_rev_table[seq_comp_table[bam1_seqi(seq, i)]];
	}
    }
    else{
	for (i = 0; i < len; ++i) {
	    buf[i] = bam_nt16_rev_table[bam1_seqi(seq, i)];
	}
    }

    sv_setpvn(scalar, buf, len);
    free(buf);
}

char *cigar_nt8_rev_table = "MIDNSHP=X";

int convert_cigar(SV* ref, int len) {
    SV* scalar = SvRV(ref);
    uint32_t* cigar = (uint32_t*)SvPV_nolen(scalar);
    
    //each cigar entry can be up to 10 characters long
    char* buf = (char*)malloc((10*len+1)*sizeof(char));
    
    //process each cigar entry
    int i;
    int op_len;
    char* cur = buf;
    char const* end = buf + sizeof(buf);
    for (i = 0; i < len; ++i) {
	//add op_len string (op_len = cigar[i] >> 4);
	cur += snprintf(cur, end - cur, "%d", cigar[i] >> 4);

	//add op to string
	*cur = cigar_nt8_rev_table[cigar[i] & 15];
	cur++;
    }

    sv_setpvn(scalar, buf, cur - buf);
    free(buf);
}

int fix_qual(SV* ref, int len, int reverse, int q64) {    
    //note these qualities are already expected to be phred+33
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    //char* buf = (char*)malloc((len+1)*sizeof(char));
    //buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(q64 == 2){ //recalculate solexa encoding
	int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 59) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-64;
	    qual[i] = 33 + qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
	}
    }
    else if(q64 != 0){ //other illumina encoding
	for (i = 0; i < len; ++i){
	    if(qual[i] < 64) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] -= 31;
	}
    }

    if (reverse != 0) {
	uint8_t b;
	for (i = 0; i < len/2; ++i){
	    if(qual[i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    if(qual[len-1-i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");

	    //swap
	    b = qual[i];
	    qual[i] = qual[len-1-i];
	    qual[len-1-i] = b;
	}
	sv_setpvn(scalar, (char*)qual, len);
    }
    else{
	for (i = 0; i < len; ++i) {
	    if(qual[i] > 75) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	}
    }
}

int convert_qual(SV* ref, int len, int reverse, int q64) {    
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    char* buf = (char*)malloc((len+1)*sizeof(char));
    buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(q64 == 2){ //recalculate solexa encoding
	int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 26) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-31;
	    qual[i] = qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
	}
    }
    else if(q64 != 0){ //other illumina encoding
	for (i = 0; i < len; ++i){
	    if(qual[i] < 31) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] -= 31;
	}
    }

    if (reverse != 0) {
	for (i = 0; i < len; ++i){
	    if(qual[i] > 42) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    buf[len - i - 1] = 33 + qual[i];
	}
    }
    else{
	for (i = 0; i < len; ++i) {
	    if(qual[i] > 42) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	    buf[i] = 33 + qual[i];
	}
    }

    sv_setpvn(scalar, buf, len);
    free(buf);
}

short reg2bin(int beg, int end){
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
    return 0;
}
