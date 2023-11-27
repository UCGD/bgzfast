#!/usr/bin/perl
use forks;
use forks::shared;

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil :sys_wait_h);
use Compress::Zlib;
use Cwd qw(abs_path);
use File::Copy;
use File::Spec;
use File::Which;
use File::Basename;
use File::Temp qw(tempfile);
use Fcntl qw(:DEFAULT);
use IPC::Shareable qw(:lock);
use Storable qw(nfreeze thaw);
use List::Util qw(shuffle);
use Perl::Destruct::Level level => 2; #helps with memory management

#constants for compression header
use constant KEY1 => pack('H8', '1f8b0804'); #4 byte static header key for BGZIP format
use constant KEY2 => pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes) for BGZIP format
use constant EOF_BLOCK => pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000'); #BGZF EOF

our @CLEANUP; #data to clean up
our $LFS;
BEGIN {
    binmode(STDIN);
    binmode(STDOUT);
    select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
    select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

    $SIG{INT} = sub {exit(130)};

    mkdir("$ENV{HOME}/.Inline") if(! -d "$ENV{HOME}/.Inline");

    $LFS = File::Which::which('lfs');
}

#use Inline (Config => DIRECTORY => File::Spec->tmpdir());
use Inline qw(C);

my ($exe) = $0 =~ /([^\/]+)$/;
my $usage = "
Usage:

     $exe <bam_file>

     Converts a BAM file to FASTQ format. Files sorted on query name convert
     faster but are not required. Paired end reads with missing pairs will be
     ignored but produce a warning. You can optionally save these missing pairs
     using the -fq3 option. Single end reads mixed in with paired end reads
     produces a fatal error unless you provide a location to store single end
     reads using the -fq3 option. Setting -fq3 to /dev/null will both throw
     aways unpaired reads as well as supress warnings and errors related to
     single end reads and missing pairs.

Options:

     fq1         <PATH>   Output FASTQ for first pair (both pairs on interleave)

     fq2         <PATH>   Output FASTQ for second pair (disabled on interleave)

     fq3         <PATH>   Output FASTQ for unmatched pairs and single-end reads

     interleave           Write output as interleaved fastq (disables -fq2)

     rg          <DIR>    Auto-generates values for -fq1, -fq2, and -fq3 using the
                          BAM read group info. If multiple read groups exist, a
                          separate -fq1, -fq2, and -fq3 set will be generated for
                          each read group. If -rg is given a direcory, then files
                          will be generated in that path. Default is to current
                          working directory.

     align_file  <PATH>   Write a FastQForward align info file to the given path

     gzip|z               Compress fastq output files using BGZF (gzip readable)

     validate|v           In memory decompression and validatation of all output
                          files post-FASTQ conversion.

     fix                  Auto-detect if q64, solexa, or restore should be used.

     q64                  Fix Phred+64 BAM qaulity scores (should be Phred+33)

     solexa               Fix Solexa odds based qaulity scores (should be Phred+33)

     restore     <TAG>    Restore original quality string using values from the
                          given BAM tag (Example: -restore OQ)

     cpus|c      <INT>    CPUs to use for conversion

     stripes     <INT>    Lustre stripes to use. Only valid for Lustre file
                          systems. Defaults to 4.

     help|?               Prints this usage statement

";

#get options from command line
my @argv = @ARGV; #backup
my %O;
GetOptions(\%O,
	   "cpus|c=i",
	   "fq1|fq=s",
	   "fq2=s",
	   "fq3=s",
	   "gzip|z",
	   "fix",
	   "q64",
	   "solexa",
	   "restore=s",
	   "validate|v",
	   "id|rg:s",
	   "stripes=i",
	   "align_file:s",
	   "interleave",
	   "test",
           "help|?") or die("ERROR: Problem with command line arguments\n");
$O{stripes} ||= 4;
$O{q64} ||= 0;
$O{q64} = 2 if($O{solexa}); #solexa goes in the $q64 value
$O{cpus} = 1 if(!$O{cpus} || $O{cpus} < 1);
delete($O{fq2}) if($O{interleave});

#help message
if($O{help}){
    print $usage;
    exit(0);
}

#make sure the id option didn't eat the file name
if($O{id} && !@ARGV && $O{id} =~ /\.bam$/){
    push(@ARGV, $O{id});
    $O{id} = '';
}

#make sure the align_file option didn't eat the file name
if($O{align_file} && !@ARGV && $O{align_file} =~ /\.bam$/){
    push(@ARGV, $O{align_file});
    $O{align_file} = '';
}

#get file
my $file = $ARGV[0];
if(!$file){
    print $usage;
    die "ERROR: Failure to provide input file\n";
}
die "ERROR: The file $file does not exist\n" if(! -f $file);
die "ERROR: Cannot specify fq2 without fq\n" if($O{fq2} && !$O{fq1});
my ($filename) = $file =~ /([^\/]+)$/;
$filename =~ s/\.bam$//i;
$O{infile} = $file;
$O{infilename} = $filename;

#generate needed file names
my %outfiles;
$O{outfiles} = \%outfiles;
prep_info($file, \%O);

#make sections of right size
my $size = (stat($file))[7];
my $count = ceil($size/33554432); #32Mb sections
$count = 5000 if($count > 5000);
my $window = ceil($size/$count);
my @temp_files : shared; #holds intermediate data from threads
my @sections = shuffle(1..$count);
share(@sections); #share is non destructive for forks::shared
$O{cpus} = $count if($count < $O{cpus});

#associate files and shared variables
our @IPC_SHARE;
our $IPCHANDLE = tie(@IPC_SHARE, 'IPC::Shareable', undef, { destroy => 1, size => 1048576 });
$IPC_SHARE[0] = 0; #output buffer offset
our %IPC4FILE;
$IPC4FILE{PIPE} = $#IPC_SHARE;

#prepare output file(s)
my $osize = 4*$size;
$osize = int($osize/3) if($O{gzip});
$O{osize} = $osize;
foreach my $set (values %outfiles){
    for(my $i = 0; $i < 3; $i++){
	my $out = $set->[$i];
	next if(!$out);

	if(! -c $out){
	    $out .= '.gz' if($O{gzip} && $out !~ /\.gz$/); #add extension
	    $out .= '~'; #make temporary file extension 
	    my ($name, $dir) = fileparse($out);
	    $dir = abs_path($dir) || abs_path('./');
	    ($out = "$dir/$name") =~ s/\/+/\//g; #cleanup
	    $set->[$i] = $out; #replace

	    #optimize for lustre file system stripes
	    if($LFS && (my $max = `$LFS osts $dir 2> /dev/null | wc -l`) > 0){
		chomp($max);
		unlink($out);
		my $scount = ($O{stripes} < $max) ? $O{stripes} : $max-1; #max includes header
		system("$LFS setstripe -c $scount -S 32m $out");
	    }

	    get_handle('>', $out); #initialize
	    close_handles(); #clear handle
	    truncate($out, $osize); #estimate that outfile needs to be 4x infile
	}

	#make ipc share for each file
	push(@IPC_SHARE, 0);
	$IPC4FILE{$out} = $#IPC_SHARE;
    }
}

#run all threads on data
for(my $i = 0; $i < $O{cpus}; $i++){
    threads->new({'context' => 'scalar'}, \&thread_run,
		 $file, $size, $count, $window, \@sections, \%O, \@temp_files, $i);
}

#clean up threads
$_->join foreach(threads->list());
unlink(@temp_files);

#finish files
foreach my $set (values %outfiles){
    for(my $i = 0; $i < 3; $i++){
	my $out = $set->[$i];
	next if(!$out);
	my $ipc_share = $IPC_SHARE[$IPC4FILE{$out}]; #size of all output
	if(!$ipc_share){
	    close_handles($out);
	    $set->[$i] = undef; #make empty
	    unlink($out) if(-f $out);
	    next;
	}
	else{
	    if($O{gzip}){
		print_async([\ (eof_block())], [$out], 0);
		$ipc_share += 28; #local copy of size must be increased
	    }
	    close_handles($out);
	    truncate($out, $ipc_share) if(-f $out);
	}

	#null padding validation
	if(-f $out){
	    sysopen(my $IN, $out, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $out: $!";
	    my $lsize = (stat($out))[7]; #logical size
	    my $hole = sysseek($IN, 0, 4); #4 is SEEK_HOLE
	    close($IN);
	    die "ERROR: File has IO related corruption: $out\n" if($lsize != $hole);
	}
    }
}

#validate that bgzip files are not corrupt
if($O{gzip} && $O{validate}){
    foreach my $set (values %outfiles){
	for(my $i = 0; $i < 3; $i++){
	    my $out = $set->[$i];
	    next if(!$out || ! -f $out);

	    #make sections of right size
	    my $size = (stat($out))[7];
	    my $count = ceil($size/33554432); #32Mb sections
	    my $window = ceil($size/$count);
	    my @sections = (1..$count);
	    share(@sections); #share is non destructive for forks::shared

	    #validate end of file
	    my $IN = get_handle('<', $out);
	    seek($IN, -28, 2); #last 28 bytes
	    my $last;
	    my $stat = read($IN, $last, 28, 0);
	    die "ERROR: Output file is missing the BGZF EOF marker: $out\n" if($stat != 28 || $last ne eof_block());
	    close_handles();

	    #run all threads on data
	    for(my $i = 0; $i < $O{cpus}; $i++){
		threads->new({'context' => 'scalar'}, \&thread_validate,
			     $file, $size, $count, $window, \@sections, \%O, $i);
	    }

	    #clean up threads
	    $_->join foreach(threads->list());	
	}
    }
    close_handles(); #clear handles
}

#move file to final location and build alignment file
foreach my $set (values %outfiles){
    for(my $i = 0; $i < 3; $i++){
	my $out = $set->[$i];
	next if(!$out || ! -f $out);
	
	(my $new = $out) =~ s/\~$//;
	File::Copy::move($out, $new) if($new ne $out);
	$set->[$i] = $new; #replace
	
	#add to alignment file entry
	next unless($O{align_file});
	my $OUT = get_handle('>', $O{align_file}) if($O{align_file});

	(my $tag = $set->[3]) =~ s/\t/\\t/g; #rg safe text tag
	if($i == 0 && $O{interleave}){
	    print $OUT "Files=$set->[0];Type=INTERLEAVED;RG=$tag\n"; 
	}
	elsif($i == 1 && !$O{interleave}){
	    print $OUT "Files=$set->[0],$set->[1];Type=PAIRED;RG=$tag\n"; 
	}
	elsif($i == 2){
	    print $OUT "Files=$set->[2];Type=UNPAIRED;RG=$tag\n";
	}
    }
}
close_handles(); #clear handles

#finished
exit(0);

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------

sub prep_info {
    my $file = shift;
    my $O = shift;

    #get read group from header
    my $header = get_header($file);
    my @rgs = grep {/^\@RG\t/} split(/\n/, $header->{text});
    
    #no read groups found
    if(!@rgs){
        warn "WARNING: No read groups specified in BAM.\n".
	    "Values will be auto-generated based off the input file name\n";
    }
    my $filename = $O->{infilename};
    my ($fullname) = $file =~ /([^\/]+)$/;
    push(@rgs, "\@RG\tID:$filename"); #always add read group for unknown
    @rgs = do {my %seen; grep { !$seen{$_}++ } @rgs}; #uniq

    #parse and split values
    my $sample;
    my %hex_seen;
    for(my $i = 0; $i < @rgs; $i++){
	my $rg = $rgs[$i];
	my %rg = map {/^([^\:]+)\:(.*)$/} grep {!/^\@RG/} split(/\t/, $rg);
	$rg{_hex} = unpack('H4', pack('S<', crc32($rg))); #helps force uniqueness	
	$rg{_old} = $rg;

	#fix repeated hexstrings
	while($hex_seen{$rg{_hex}}++){
	    my $crc = unpack('S<', pack('H4', $rg{_hex})); #convert back to int
	    $rg{_hex} = unpack('H4', pack('S<', $crc + 1)); #iterate and convert to hex
	}

	#WashU data
	if(defined($rg{SM}) && defined($rg{LB}) && $rg{LB} =~ /^\"?$rg{SM}\-.*lib\d+\"?$/){
	    my ($washu_prefix) = $rg{SM} =~ /^(H_[A-Z]{2}\-)[^\-]+\-[^\-]+$/; #version1 prefix
	    ($washu_prefix) = $rg{SM} =~ /^(H_[A-Z]{2}\-)[^\-]+$/ if(!$washu_prefix); #version2 prefix
	    if($washu_prefix){
		$washu_prefix = quotemeta($washu_prefix);
		$rg{SM} =~ s/^$washu_prefix//;
		($rg{_lane} = $rg{PU}) =~ s/^.*\.(\d+)$/L$1/;
	    }
	}
	$rgs[$i] = \%rg; #replace
	$sample = $rg{SM} if(!defined($sample) && defined($rg{SM}));
    }
    $sample = $filename if(!defined($sample)); #assumes all data in file is from the same sample

    #fix values
    for(my $i = 0; $i < @rgs; $i++){
        my %rg = %{$rgs[$i]};

	#fill in missing values
	$rg{ID} = $rg{_hex} if(!defined($rg{ID}));
	$rg{SM} = $sample if(!defined($rg{SM})); #assumes all data is from the same sample
	$rg{LB} = "$rg{SM}-lib1" if(!defined($rg{LB}));
	$rg{PL} = 'ILLUMINA' if(!defined($rg{PL}));
	$rg{PL} = uc($rg{PL}); #must be all caps
	$rg{PU} = $rg{PL}."-1" if(!defined($rg{PU}));
	($rg{_lane} = $rg{PU}) =~ s/^${rg{PL}}[\-_]L?(\d+)/L$1/i if(!defined($rg{_lane}));
	
	#make final tag
	$rg{_rg} = "\@RG\tID:$rg{ID}\tSM:$rg{SM}\tLB:$rg{LB}\tPL:$rg{PL}\tPU:$rg{PU}";
	$rgs[$i] = \%rg; #replace
    }

    #evaluate alignments
    if(defined($O->{fix})){
	my %stats;
	my $size = (stat($file))[7];
	my $block = $header->{block}; #piece of alignment stuck at end of header
	my $pos = $header->{offset};
	my $boff = 0;
	my $IN = get_handle('<', $file);
	seek($IN, $pos, 0);
	while($pos < $size || length($block) < $boff){
	    inflate(readblock($IN), \$block);
	    $pos = tell($IN);
	    $boff = eval_alignments(\$block, $boff, \%stats);
	    last if($stats{count} >= 50000); #only evaluate 50000 reads
	}
	my $restore = ($stats{tag}{OQ}) ? 'OQ' : '';
	my $q64 = qual_detectformat($stats{min}, $stats{max}, $stats{mean});
	$O->{restore} = $restore if(!$O->{restore});
	$O->{q64} = $q64 if(!$O->{q64});
    }

    #make alignment info file
    if(defined($O->{align_file}) && $O->{align_file} eq ''){
	$O->{align_file} = "./$fullname.bam2fastq.info";
	print STDERR "##NOTE: Alignment info will be sent to $O->{align_file}\n";
    }        

    #make output files based on ID
    if(defined($O->{id})){ #files for each read group
	my $dir = $O->{id};
	$dir = "./$fullname.bam2fastq" if($dir eq '');
	$dir =~ s/\/+$//; #remove trailing '/'
	mkdir($dir) if(!-d $dir);
	$O->{id} = $dir;

	warn "WARNING: Multiple read groups specified in BAM.\n".
	    "Separate output file will be created for each.\n" if(@rgs > 2);
	
	for(my $i = 0; $i < @rgs; $i++){
	    my $rg = $rgs[$i];
	    my ($rg_id, $hex, $lane, $tag) = ($rg->{SM}, $rg->{_hex}, $rg->{_lane}, $rg->{_rg});

	    #make file names
	    my $name = "$rg_id\_$hex\_$lane";
	    my $out1 = "$dir/$name\_1.fq";
	    my $out2 = "$dir/$name\_2.fq";
	    my $out3 = "$dir/$name.fq";
	    if($O->{interleave}){
		$out1 = "$dir/$name\_interleaved.fq";
		$out2 = undef;
	    }

	    #add extension
	    if($O->{gzip}){
		foreach my $out ($out1, $out2, $out3){
		    next if(!length($out));
		    $out .= '.gz' if($out !~ /\.gz$/);
		}
	    }

	    #let user know where the output will be
	    if($i != $#rgs){
		print STDERR "##NOTE: RG $rg->{ID} will be sent to $out1";
		print STDERR ", $out2," if($out2);
		print STDERR " and $out3\n";
	    }
	    else{
		print STDERR "##NOTE: Reads with no RG tag will be sent to $out1";
		print STDERR ", $out2," if($out2);
		print STDERR " and $out3\n";
	    }
	    $O->{outfiles}{$rg->{ID}} = [$out1, $out2, $out3, $tag];
	}

	#clear
	delete($O->{fq1});
	delete($O->{fq2});
	delete($O->{fq3});
	
	exit(0) if($O->{test});
    }
    else{ #just use first readgroup tag
	my ($out1, $out2, $out3) = ($O->{fq1}, $O->{fq2}, $O->{fq3});
	$out2 = undef if($O->{interleave});
	if($O->{gzip}){
	    foreach my $out ($out1, $out2, $out3){
		next if(!length($out));
		$out .= '.gz' if($out !~ /\.gz$/);
	    }
	}    
	$O->{outfiles}{$filename} = [$out1, $out2, $out3, $rgs[0]{_rg}]; 

	print STDERR "##NOTE: Only the first RG tag will be used for the align info file,\n".
	    "and it will apply to all reads reguardless of the original read group\n"
	    if(defined($O->{align_file}));
    }
}

#get header from bam
sub get_header {
    my $file = shift;

    #declare variables
    my $header;
    my $header_off;
    my $header_block;
    
    #open input file
    my $IN = get_handle('<', $file);
    
    #get header
    inflate(readblock($IN), \$header_block);
    
    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
	if($magic ne "BAM\1");
    
    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
	inflate(readblock($IN), \$header_block);
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
	    inflate(readblock($IN), \$header_block);
	}
	my $l_name = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	
	while(length($header_block)-$header_off < $l_name+4){ #grow if needed
	    inflate(readblock($IN), \$header_block);
	}
	my $name = substr($header_block, $header_off, $l_name);
	$header_off += $l_name;
	$name =~ s/\0$//g; #remove null padding
	$header->{ref}[$i]{name} = $name;
	
	my $l_ref = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	$header->{ref}[$i]{l_ref} = $l_ref;
    }
    substr($header_block,0,$header_off,''); #chop off header from block
    $header_off = tell($IN); #make offset be position of first block after header

    $header->{offset} = $header_off; #offset of first block imediately following the header
    $header->{block}  = $header_block; #piece of alignment accidentally stuck inside header block

    return $header;
}

sub eval_alignments {
    my $data  = shift; #ref
    my $block_offset = shift || 0;
    my $stats = shift;

    #collect and process bam alignments
    my $data_len = length($$data);
    while($block_offset < $data_len - 4){ #continue until near end of block
        my $substr_off = $block_offset;

        my ($block_size) = unpack('l<', substr($$data, $substr_off, 4));
        $substr_off += 4;
        last if($substr_off+$block_size > $data_len); #current alignment is spit across next block
        $block_offset = $substr_off+$block_size; #end of current alignmnet

        #get valuse
        #($refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen)
        #    = unpack('l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 32));

        #flag_nc processing into sub values
        my ($flag_nc) = unpack('L<', substr($$data, $substr_off+12, 4));
        my $flag = ($flag_nc>>16); #flag is probably not usful here
        next if($flag & 2816); #skip secondary, supplemental, and vendor failed alignments
        my $n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1

        #get sequence length
        my ($l_seq) = unpack('l<', substr($$data, $substr_off+16, 4));

        #get bin_mq_nl and process into sub values
        my ($bin_mq_nl) = unpack('L<', substr($$data, $substr_off+8, 4));
        #my $bin = ($bin_mq_nl>>16);
        #my $mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1
        my $l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
        $substr_off += 32;
        $substr_off += $l_read_name; #jump past read name

        #move past cigar string
        $substr_off += $n_cigar_op*4;

        #jump past seq string
        my $count = int(($l_seq+1)/2);
        $substr_off += $count;

        #get quality string
        my $qual = substr($$data, $substr_off, $l_seq);
        $substr_off += $l_seq;

        #restore quality value from tags
	my $oq;
        while($substr_off < $block_offset){
            my $tag =  substr($$data, $substr_off, 2);
            $stats->{tag}{$tag}++;
            $substr_off += 2;
            my $val_type = substr($$data, $substr_off, 1);
            $substr_off += 1;

            my $val_len = 0;
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
                my $sub_type = substr($$data, $substr_off, 1);
                $substr_off += 1;
                my $sub_len = unpack('l<', substr($$data, $substr_off, 4));
                $substr_off += 4;

                if   ($sub_type eq 'c'){ $val_len = 1*$sub_len }
                elsif($sub_type eq 'C'){ $val_len = 1*$sub_len }
                elsif($sub_type eq 's'){ $val_len = 2*$sub_len }
                elsif($sub_type eq 'S'){ $val_len = 2*$sub_len }
                elsif($sub_type eq 'i'){ $val_len = 4*$sub_len }
                elsif($sub_type eq 'I'){ $val_len = 4*$sub_len }
                elsif($sub_type eq 'f'){ $val_len = 4*$sub_len }
            }

            if($tag eq 'OQ'){
                die "ERROR: Wrong datatype for original qaulity value\n" if($val_type ne 'Z');
                die "ERROR: Restored quality value does not match sequence length\n" if($val_len-1 != $l_seq);
                my $value = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
                $substr_off += $val_len;
		$oq = 1;

                #replace qual
                $qual = $value;
            }
            else{ #skip past value since it is not the right one
                $substr_off += $val_len;
            }
        }

	#convert Phred+0 to Phred+33
	if(!$oq){
	    qual_0_to_33(\$qual, $l_seq, 0, 0);
	}

        #evaluate qual range
        qual_stats(\$qual, $l_seq, $stats);
    }

    return $block_offset;
}

sub qual_detectformat {
    my $min = shift;
    my $max = shift;
    my $mean = shift;

    if(!defined($min) || !defined($max)){
        die "ERROR: Both min and max must be provided to determine format\n";
    }
    elsif(33 <= $min && $max <= 78){ #Phred+33 (allow qual 0 to 45 plus 33)
        return 0; #'Phred+33';
    }
    elsif(64 <= $min && $max <= 109){ #Phred+64 (allow qual 0 to 45 plus 64)
        return 1; #'Phred+64';
    }
    elsif(59 <= $min && $max <= 109){ #Solexa+64 (allow qual -5 to 45 plus 64)
        return 2; #'Solexa+64';
    }
    #elsif(0 <= $min && $max <= 45){ #Phred+0 (allow qual 0 to 45 plus 0)
    #    return 'Phred+0';
    #}
    #elsif(31 <= $min && $max <= 76){ #Phred+64 (allow qual 0 to 45 plus 31)
    #   return 'Phred+31';
    #}
    #elsif(26 <= $min && $max <= 76){ #Solexa+64 (allow qual -5 to 45 plus 31)
    #   return 'Solexa+31';
    #}
    else{ #make best guess
        if(abs(104-$mean) < abs(73-$mean)){
            return 1; #'Phred+64';
        }
        else{
            return 0; #'Phred+33';
        }
    }
}

#convenience method to see if I am in a thread or not
{my $tid; BEGIN{$tid = threads->tid if exists $INC{'threads.pm'}};
sub is_thread {
    my $is_thread = 1 if(defined($tid) && $tid != threads->tid);
    return $is_thread;
}}

#main fastq conversion done here by each process
sub thread_run {
    my $file       = shift;
    my $size       = shift;
    my $count      = shift;
    my $window     = shift;
    my $sections   = shift;
    my $O          = shift;
    my $temp_files = shift;
    my $id         = shift;

    #fix thread signalling
    if(is_thread){
	select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
	select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

	$SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
	$SIG{INT} = sub {exit(130)};
    }

    #load parameters
    my $outfiles   = $O->{outfiles};
    my $q64        = $O->{q64};
    my $restore    = $O->{restore};
    my $gzip       = $O->{gzip};
    my $interleave = $O->{interleave};
    my $osize      = $O->{osize};
    my $cpus       = $O->{cpus};

    #declare variables
    my @forks;
    my $header;
    my $header_off;
    my $header_block;
    my $paired;
    my %buffers;
    foreach my $rg (keys %outfiles){
	my $fq1 = '';
	my $fq2 = '';
	my $fq3 = '';
	my $gz1 = '';
	my $gz2 = '';
	my $gz3 = '';

	#configure output buffers
	my $fq_buffer;
	if($interleave){
	    if($gzip){
		$fq_buffer = [\$fq1, \$fq1, \$fq3, \$gz1, \$gz1, \$gz3];
	    }
	    else{
		$fq_buffer = [\$fq1, \$fq1, \$fq3, \$fq1, \$fq1, \$fq3];
	    }
	}
	else{
	    if($gzip){
		$fq_buffer = [\$fq1, \$fq2, \$fq3, \$gz1, \$gz2, \$gz3];
	    }
	    else{
		$fq_buffer = [\$fq1, \$fq2, \$fq3, \$fq1, \$fq2, \$fq3];
	    }
	}

	#configure read buffer
	my $max = ($cpus > 100) ? $cpus : 100;
	my $r_buffer = [map { {} } (1..$max)]; #paired read buffer (array of hash refs)

	#add to buffer set
	$buffers{$rg} = [$fq_buffer, $r_buffer];
    }

    #open input file
    my $IN = get_handle('<', $file);
    
    #get header
    inflate(readblock($IN), \$header_block);
    
    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
	if($magic ne "BAM\1");
    
    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
	inflate(readblock($IN), \$header_block);
    }
    my $text = substr($header_block, $header_off, $l_text);
    $header_off += $l_text;
    $text =~ s/\0$//g; #remove null padding
    $header->{text} = $text;
    
    my $n_ref = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    $header->{n_ref} = $n_ref;
    
    #get detail on each reference
    for(my $i = 0; $i < $n_ref; $i++){
	while(length($header_block)-$header_off < 4){ #grow if needed
	    inflate(readblock($IN), \$header_block);
	}
	my $l_name = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	
	while(length($header_block)-$header_off < $l_name+4){ #grow if needed
	    inflate(readblock($IN), \$header_block);
	}
	my $name = substr($header_block, $header_off, $l_name);
	$header_off += $l_name;
	$name =~ s/\0$//g; #remove null padding
	$header->{ref}[$i]{name} = $name;
	
	my $l_ref = unpack('l<', substr($header_block, $header_off, 4));
	$header_off += 4;
	$header->{ref}[$i]{l_ref} = $l_ref;
    }
    substr($header_block,0,$header_off,''); #chop off header from block
    $header_off = tell($IN); #make offset be position of first block after header

    #storage file to hold data to avoid high RAM usage
    my ($tfh, $tfile) = tempfile("pairbam\_$id\_XXXXX", CLEANUP => 1, DIR => File::Spec->tmpdir);
    push(@CLEANUP, $tfile);

    #now read alignment sections
    while(my $section = shift @$sections){
	my $start = ($section-1) * $window;
	my $end   = $section * $window;
	$start = $header_off if($start < $header_off);
	$end = $size if($end > $size); #don't go after end

	#adjust to start of next block
	seek_next_block($IN, $start, 0);
	
	#get first block in section
	my $block = ($section == 1) ? $header_block : ''; #piece of alignment stuck at end of header
	inflate(readblock($IN), \$block);

	#find first alignment
	my $offset = tell_next_alignment(\$block, $header);
	while(! defined($offset) && tell($IN) < $start+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$block);
	    $offset = tell_next_alignment(\$block, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($block));
	substr($block, 0, $offset) = ''; #chop off leading data

	#get all blocks in section
	while(tell($IN) < $end){
	    inflate(readblock($IN), \$block);
	    bam2fastq(\$block, \%buffers, $O);
	}

	#get last block and adjust split terminating line
	my $tail = '';
	inflate(readblock($IN), \$tail);
	$offset = tell_next_alignment(\$tail, $header);
	while(! defined($offset) && tell($IN) < $end+65536 && tell($IN) < $size){ #get past short and empty blocks
	    inflate(readblock($IN), \$tail);
	    $offset = tell_next_alignment(\$tail, $header);
	}
	die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($tail));
	$block .= substr($tail, 0, $offset) if($offset); #chop off trailing data

	#convert to fastq
	bam2fastq(\$block, \%buffers, $O);

	#print section to files for each read group
	foreach my $rg (keys %buffers){
	    my $fq_buffer = $buffers{$rg}[0];
	    my ($out1, $out2, $out3, $tag) = @{$outfiles->{$rg}};
	    
	    #compress remainder in section
	    if($gzip){
		#max ISIZE is 65536 - 256
		deflate(\ (substr(${$fq_buffer->[0]}, 0, 65280)), $fq_buffer->[3])
		    while(length(${$fq_buffer->[0]}));
		deflate(\ (substr(${$fq_buffer->[1]}, 0, 65280)), $fq_buffer->[4])
		    while(length(${$fq_buffer->[1]}));
		deflate(\ (substr(${$fq_buffer->[2]}, 0, 65280)), $fq_buffer->[5])
		    while(length(${$fq_buffer->[2]}));
	    }
	    
	    #minimal validation
	    if(!$out3){
		$paired |= 1 if(length(${$fq_buffer->[5]})); #single end
		$paired |= 2 if(length(${$fq_buffer->[3]}) || length(${$fq_buffer->[4]})); #paired end
		if($paired & 1 && ($out2 || $interleave)){
		    die "ERROR: The input BAM contains single end reads. You indicated\n".
			"that you expect paired end reads, but you failed to specify a\n".
			"file for single end reads using -fq3.\n";
		}
		if($paired & 1 && $paired & 2){
		    die "ERROR: The input BAM contains both paired and single end reads,\n".
			"but you failed to specify a file for single end reads using -fq3\n";
		}
	    }
		
	    
	    print_async([$fq_buffer->[3], $fq_buffer->[4], $fq_buffer->[5]],
			[$out1, $out2, $out3], 0);
	}

	#stores buffer to avoid high memory usage
	foreach my $rg (sort keys %buffers){
	    my $r_buffer = $buffers{$rg}[1];
	    next if(!@$r_buffer);
	    
	    for(my $i = 0; $i < @$r_buffer; $i++){
		next if(! keys %{$r_buffer->[$i]}); #skip empty hashes
		my $dref = zip(nfreeze($r_buffer->[$i]));
		undef(%{$r_buffer->[$i]}); #clear memory
		print $tfh pack('S<S<Q<', $i, length($rg), length($$dref)); #(fasthashid, rg_len, buf_len)
		print $tfh $rg;
		print $tfh $$dref;
	    }
	}
    }

    #close storage file and share among threads
    close($tfh);
    push(@$temp_files, $tfile);
    sleep 1 while(@$temp_files < $cpus);

    #gather and try and to colapse remaining temp_file buffers
    my $max = ($cpus > 100) ? $cpus : 100;
    my @search = grep {$_ % $cpus == $id} (0..$max);
    foreach my $sid (@search){
	foreach my $file (@$temp_files){
	    open($tfh, "< $file");
	    my $meta;
	    while(read($tfh, $meta, 12)){
		my ($fhid, $rglen, $len) = unpack('S<S<Q<', $meta); #is this the buffer section I want
		if($fhid != $sid){ #skip past unwanted buffer sections
		    seek($tfh, $rglen+$len, 1);
		    next;
		}

		#retrieve buffer section
		my $rg;
		my $hash2;
		read($tfh, $rg, $rglen);
		read($tfh, $hash2, $len);
		$hash2 = thaw(${unzip(\$hash2)});
		my $fq_buffer = $buffers{$rg}[0];
		my $r_buffer  = $buffers{$rg}[1];
		if(!$r_buffer->[$sid]){
		    $r_buffer->[$sid] = $hash2;
		    next;
		}
		$r_buffer->[$sid] ||= {};
		my $hash1 = $r_buffer->[$sid];
		
		#find reads between hashes
		foreach my $read_name (keys %$hash2){
		    if(!exists($hash1->{$read_name})){ #copy over to local hash
			$hash1->{$read_name} = delete($hash2->{$read_name});
			next;
		    }
		    
		    if(chop($hash1->{$read_name}) == 1){ #first in pair
			chop($hash2->{$read_name});
			${$fq_buffer->[0]} .= delete($hash1->{$read_name});
			${$fq_buffer->[1]} .= delete($hash2->{$read_name});
		    }
		    else{ #second in pair
			chop($hash2->{$read_name});
			${$fq_buffer->[0]} .= delete($hash2->{$read_name});
			${$fq_buffer->[1]} .= delete($hash1->{$read_name});
		    }
		    
		    #compress
		    if($O{gzip}){
			#max ISIZE is 65536 - 256
			deflate(\ (substr(${$fq_buffer->[0]}, 0, 65280)), $fq_buffer->[3])
			    while(length(${$fq_buffer->[0]}) >= 65280);
			deflate(\ (substr(${$fq_buffer->[0]}, 0, 65280)), $fq_buffer->[4])
			    while(length(${$fq_buffer->[0]}) >= 65280);
		    }
		    
		    #print
		    if(length(${$fq_buffer->[3]}) >= 33554432 || length(${$fq_buffer->[4]}) >= 33554432){
			my ($out1, $out2, $out3, $tag) = @{$outfiles{$rg}};
			print_async([$fq_buffer->[3], $fq_buffer->[4]], [$out1, $out2], 0);
		    }
		}
		undef(%{$hash2}); #coerce perl to clear RAM where possible
	    }
	    close($tfh);
	}

	#handle reads missing mates and clear remaining buffer
	foreach my $rg (keys %buffers){
	    my $bset = $buffers{$rg};
	    my $fq_buffer = $bset->[0];
	    my $r_buffer = $bset->[1];
	    
	    my $outs = $outfiles{$rg};
	    my ($out1, $out2, $out3, $tag) = @{$outs};
	    
	    $r_buffer->[$sid] ||= {}; #just in case
	    my $hash1 = $r_buffer->[$sid];
	    
	    foreach my $read_name (keys %$hash1){
		my $read = delete($hash1->{$read_name});
		chop($read);
		if($out3){
		    ${$fq_buffer->[2]} .= $read;
		    
		    #compress
		    if($O{gzip}){
			#max ISIZE is 65536 - 256
			deflate(\ (substr(${$fq_buffer->[2]}, 0, 65280)), $fq_buffer->[5])
			    while(length(${$fq_buffer->[2]}) >= 65280);
		    }
		    print_async([$fq_buffer->[5]], [$out3], 0) if(length(${$fq_buffer->[5]}) >= 33554432);
		}
		else{
		    die "ERROR: paired read $read_name did not have a mate\n".
			"No -fq3 file was specified to hold orphan mates\n";
		}
	    }

	    #compress what's left
	    if($O{gzip}){
		#max ISIZE is 65536 - 256
		deflate(\ (substr(${$fq_buffer->[0]}, 0, 65280)), $fq_buffer->[3])
		    while(length(${$fq_buffer->[0]}));
		deflate(\ (substr(${$fq_buffer->[1]}, 0, 65280)), $fq_buffer->[4])
		    while(length(${$fq_buffer->[1]}));
		deflate(\ (substr(${$fq_buffer->[2]}, 0, 65280)), $fq_buffer->[5])
		    while(length(${$fq_buffer->[2]}));
	    }
	    print_async([$fq_buffer->[3], $fq_buffer->[4], $fq_buffer->[5]],
			[$out1, $out2, $out3], 0);
	}
    }   
    
    #free resources
    undef %buffers;
    close_handles();
    undef @CLEANUP;

    return 1;
}

our %HANDLES;
sub get_handle {
    my $mode = shift;
    my $file = shift;

    if(!$HANDLES{$mode}{$file}){
	my $mode2;
	if($mode eq '<'){
	    $mode2 = O_RDONLY|O_BINARY;
	}
	elsif($mode eq '>'){
	    $mode2 = O_WRONLY|O_CREAT|O_TRUNC|O_BINARY;
	}
	elsif($mode eq '>>'){
	    $mode2 = O_WRONLY|O_CREAT|O_BINARY;
	}
	elsif($mode eq '+<'){
	    $mode2 = O_RDWR|O_CREAT|O_BINARY;
	}
	else{
	    die "ERROR: Unsupported mode";
	}

	sysopen(my $FH, $file, $mode2) or return undef;
	$HANDLES{$mode}{$file} = $FH;
    }

    return $HANDLES{$mode}{$file};
}
sub close_handles {
    if(@_){
	foreach my $file (@_){
	    foreach my $mode (keys %HANDLES){
		next unless(exists($HANDLES{$mode}{$file}));
		my $FH = delete($HANDLES{$mode}{$file});
		close($FH);
	    }
	}
    }
    else{
	close($_) foreach(map {values %$_} values %HANDLES);
	undef %HANDLES;
    }
}

sub print_async {
    my $refs     = shift;
    my $outfiles = shift;
    my $async    = shift;

    #lock and change file size if necessary
    my $ipc_handle = $IPCHANDLE;
    if($async){
	return unless($ipc_handle->shlock(LOCK_EX|LOCK_NB));
    }
    else{
	sleep 0.1 while(!$ipc_handle->shlock(LOCK_EX));
    }

    my @locs;
    for(my $i = 0; $i < @$refs; $i++){
	my $ref = $refs->[$i];
	next if(!length($$ref));

	my $out = $outfiles->[$i];
	next if(!$out);

	my $share_id = $IPC4FILE{$out};
	my $ipc_share = $IPC_SHARE[$share_id];

	#output first pair or interleaved pairs
	$locs[$i] = $ipc_share;
	$ipc_share += length($$ref);
	$IPC_SHARE[$share_id] = $ipc_share;
	if(-c $out){
	    while(length($$ref)){ #keep trying until it prints everything
		my $stat = syswrite(STDOUT, $$ref);
		$stat = syswrite(STDOUT, $$ref) if(!$stat);
		die "Writing error: $!" if(!defined($stat));
		substr($$ref, 0, $stat, ''); #reset
	    }
	}
	else{
	    truncate($out, $ipc_share) if($ipc_share > $osize); #grow
	}
    }
    $ipc_handle->shunlock;

    #write data once lock is no longer needed
    for(my $i = 0; $i < @$refs; $i++){
	my $ref = $refs->[$i];
	my $out = $outfiles->[$i];

	#output first pair or interleaved pairs
	next if(! length($$ref));

	my $loc = $locs[$i];
	die "ERROR: Cannot write unless outfile specified" if(!$out);
	
	my $OUT = get_handle('+<', $out);
	sysseek($OUT, $loc, 0);
	while(length($$ref)){ #keep trying until it prints everything
	    my $stat = syswrite($OUT, $$ref);
	    while(!defined($stat)){
		close_handles(); #just close them all
		$OUT = get_handle('+<', $out);
		if(!$OUT){
		    sleep 10;
                    next;
		}
		sysseek($OUT, $loc, 0);
		$stat = syswrite($OUT, $$ref);
		if(!defined($stat) && $! =~ /Interrupted system call/){
		    sleep 10;
		    next;
		}
		last;
	    }
	    die "Writing error: $!" if(!defined($stat));
	    $loc += $stat;
	    substr($$ref, 0, $stat, ''); #reset
	}
    }
}
    
{#predeclare values for efficienty
my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$count,$next_refID,$next_pos,$tlen,
    $bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,$seq,$qual,$rg,$oq);
my ($tag, $val_type, $val_len, $sub_type, $sub_len);
my ($cl_ref,$cl_query,$match, $read);
my ($substr_off, $data_len, $block_offset, $qual_off);
my ($data, $buffers, $O, $q64, $restore, $gzip, $filename);
my ($last, $bset, $fq_buffer, $r_buffer, $buf_count, $r_buf_bins, $id);
my ($oq_warn);

sub bam2fastq {
    $data      = shift; #ref
    $buffers   = shift; #ref
    $O         = shift; #ref
    $q64       = $O->{q64};
    $restore   = $O->{restore};
    $gzip      = $O->{gzip};
    $filename = $O->{infilename};

    #collect and process bam alignments
    $data_len = length($$data);
    $last = $filename;
    $bset = $buffers->{$filename};
    $fq_buffer  = $bset->[0];
    $r_buffer   = $bset->[1];
    $buf_count = scalar(keys %$buffers);
    $r_buf_bins = scalar(@{$r_buffer});
    $block_offset = 0;
    while($block_offset < $data_len - 4){ #continue until near end of block
	$substr_off = $block_offset;
	
	($block_size) = unpack('l<', substr($$data, $substr_off, 4));
	$substr_off += 4;
	last if($substr_off+$block_size > $data_len); #current alignment is spit across next block
	$block_offset = $substr_off+$block_size; #end of current alignmnet
	
	($refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen)
	    = unpack('l<l<L<L<l<l<l<l<', substr($$data, $substr_off, 32));
	$substr_off += 32;
	
	#bin_mq_nl processing into sub values
	#$bin = ($bin_mq_nl>>16);
	#$mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1 
	$l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	
	#flag_nc processing into sub values
	$flag = ($flag_nc>>16); #flag is probably not usful here
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	next if($flag & 2816); #skip secondary, supplemental, and vendor failed alignments
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name-1); #-1 to drop null
	$substr_off += $l_read_name;
	
	#process cigar into sub values
	if($n_cigar_op){
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
	    if($pos >= 0 && $cl_query != $l_seq){ #hard masking.  I've lost sequence
		warn "WARNING: $read_name has lost sequence because of hard masking of reads\n";
	    }
	}
	
	#process seq string
	$count = int(($l_seq+1)/2);
	$seq = substr($$data, $substr_off, $count);
	$substr_off += $count;
	seq_bam2string(\$seq, $l_seq, ($flag & 16));
	
	#process quality string
	$qual_off = $substr_off;
	$substr_off += $l_seq;

	#process attributes tags
	if($buf_count > 1 || $restore){
	    $oq = 0;
	    $rg = '';
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
		
		#get read group tag
		if($buf_count > 1 && $tag eq 'RG'){
		    die "ERROR: Wrong datatype for read group tag\n" if($val_type ne 'Z');
		    $rg = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
		    $substr_off += $val_len;
		    last if(!$restore || length($oq)); #short circuit
		}

		if($restore && $tag eq $restore){
		    die "ERROR: Wrong datatype for restored qaulity value\n" if($val_type ne 'Z');
		    die "ERROR: Restored quality value does not match sequence length\n" if($val_len-1 != $l_seq);
		    $qual = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
		    qual_33_fix(\$qual, $l_seq, $q64, ($flag & 16));
		    $substr_off += $val_len;
		    $oq = 1; #found OQ
		    last if($buf_count <= 1 || length($rg)); #short circuit loop
		}
		else{ #skip past value since it is not the right one
		    $substr_off += $val_len;
		}
	    }
	}
	
	if(!$oq){
	    if($restore && !$oq_warn){
		warn "WARNING: Original quality value not found for some reads\n";
		$oq_warn++;
	    }

	    $qual = substr($$data, $qual_off, $l_seq);
	    qual_0_to_33(\$qual, $l_seq, $q64, ($flag & 16));
	}
	
	#use read group to load correct buffers
	$rg = $filename if(!length($rg)); #empty read group
	if($rg ne $last){
	    #compress
	    if($gzip){
		#max ISIZE is 65536 - 256
		deflate(\ substr(${$fq_buffer->[0]}, 0, 65280), $fq_buffer->[3])
		    while(length(${$fq_buffer->[0]}) >= 65280);
		deflate(\ substr(${$fq_buffer->[1]}, 0, 65280), $fq_buffer->[4])
		    while(length(${$fq_buffer->[1]}) >= 65280);
		deflate(\ substr(${$fq_buffer->[2]}, 0, 65280), $fq_buffer->[5])
		    while(length(${$fq_buffer->[2]}) >= 65280);
	    }
	    
	    $last = $rg;
	    $bset = $buffers->{$rg};
	    $fq_buffer  = $bset->[0];
	    $r_buffer   = $bset->[1];
	    $rg = ''; #reset
	}
	
	#group pairs together
	if(!($flag & 1)){ #read not paired
	    ${$fq_buffer->[2]} .= "\@$read_name\n$seq\n+\n$qual\n";
	}
	else{ #pairs
	    $id = fasthash($read_name) % $r_buf_bins;
	    if(exists($r_buffer->[$id]{$read_name})){
		$read = delete($r_buffer->[$id]{$read_name});
		chop($read);
		if($flag & 64){ #first in pair
		    ${$fq_buffer->[0]} .= "\@$read_name/1\n$seq\n+\n$qual\n";
		    ${$fq_buffer->[1]} .= $read;
		}
		else{
		    ${$fq_buffer->[0]} .= $read;
		    ${$fq_buffer->[1]} .= "\@$read_name/2\n$seq\n+\n$qual\n";
		}
	    }
	    else{
		if($flag & 64){ #first in pair
		    $r_buffer->[$id]{$read_name} = "\@$read_name/1\n$seq\n+\n$qual\n1"; #pair number is last character of string
		}
		else{
		    $r_buffer->[$id]{$read_name} = "\@$read_name/2\n$seq\n+\n$qual\n2"; #pair number is last character of string
		}
	    }
	}
    }

    #compress
    if($gzip){
	#max ISIZE is 65536 - 256
	deflate(\ substr(${$fq_buffer->[0]}, 0, 65280), $fq_buffer->[3])
	    while(length(${$fq_buffer->[0]}) >= 65280);
	deflate(\ substr(${$fq_buffer->[1]}, 0, 65280), $fq_buffer->[4])
	    while(length(${$fq_buffer->[1]}) >= 65280);
	deflate(\ substr(${$fq_buffer->[2]}, 0, 65280), $fq_buffer->[5])
	    while(length(${$fq_buffer->[2]}) >= 65280);
    }

    $$data = substr($$data, $block_offset); #remove what I processed off of the block

    return;
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
    while(1){
	#match first static key
	$offset = index($data, KEY1, $offset);
	return seek($FH, $size, 0) if($offset == -1);
	
	#second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
	last if(substr($data, $offset+10, 6) eq KEY2);
	
	$offset += 4; #jump over key1 match and try again
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

    die "ERROR: bsize is too large, BGZF file may be corrupt\n" if($bsize >= 65536);
    
    #read compression block and footer
    my $c_data_size = $bsize-$xlen-19; #compression block size
    $stat = read($FH, $data, $c_data_size + 8, length($data)); #the +8 is for the footer 
    die "ERROR: Could not read compression block\n"
	if($stat != $c_data_size + 8);

    #minimal validation
    my ($isize) = unpack('L<', substr($data, -4, 4));
    die "ERROR: isize is too large, BGZF file may be corrupt\n" if($isize > 65536);

    return \$data;
}

sub inflate {
    my $in_buffer  = $_[0]; #reference
    my $out_buffer = $_[1]; #reference

    return if(!$in_buffer);
    return if(!length($$in_buffer));

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 131072);
    die "ERROR: Failed to create zlib inflation object with status: $stat" if($stat);

    substr($$in_buffer, 0, 18, ''); #drop header without validation
    $$out_buffer .= scalar($i_obj->inflate($in_buffer)); #destructive
    die "ERROR: Trailing garbage in compression block\n" if(length($$in_buffer) != 8);

    #validate the length and crc32
    my ($crc32, $isize) = unpack('L<L<', $$in_buffer); #footer
    if($isize != $i_obj->total_out()){ #size does not match
        die "ERROR: The expected ISIZE of the uncompressed block does not match\n";
    }
    if($crc32 != crc32(substr($$out_buffer, -$isize, $isize))){ #crc32 does not match
        die "ERROR: The expected CRC32 of the uncompressed block does not match\n";
    }
    $$in_buffer = ''; #empty

    return;
}

sub deflate {
    my $in_buffer = $_[0]; #reference
    my $out_buffer = $_[1]; #reference
    my $level = $_[2];

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 131072, -Level => 4);
    die "ERROR: Failed to create zlib deflation object with status: $stat\n" if($stat);

    my $offset = length($$out_buffer);
    $$out_buffer .= pack('H36', '1f8b08040000000000ff0600424302000000'); #header (BSIZE=0)
    $$out_buffer .= scalar($d_obj->deflate($in_buffer)).scalar($d_obj->flush); #compressed data
    $$out_buffer .= pack('V', crc32($in_buffer)); #CRC32
    $$out_buffer .= pack('V', $d_obj->total_in()); #ISIZE
    substr($$out_buffer, $offset+16, 2, pack('v', $d_obj->total_out+6+19)); #set final BSIZE
    $$in_buffer = ''; #destroy input buffer

    return;
}

sub zip {
    my $iref  = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $oref  = $_[1];

    if(!$oref){
        my $string = '';
        $oref = \$string;
    }
    elsif(!ref($oref)){
        $oref = \ ($_[1]);
    }
    $$oref = '' if(!defined($$oref));

    return $oref if(!length($$iref));

    #compress data from string to string
    while(length($$iref)){
        deflate(\ (substr($$iref, 0, 65280)), $oref);
    }

    return $oref;
}

sub unzip {
    my $iref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $oref = $_[1];

    if(!$oref){
        my $string = '';
        $oref = \$string;
    }
    elsif(!ref($oref)){
        $oref = \ ($_[1]);
    }

    return $oref if(!length($$iref));

    #decompress data from string to string
    my $bsize; #predeclare
    while(length($$iref) >= 28){
        #validate BGZF header
        die "ERROR: Does not appear to be a BGZF file\n"
            if(substr($$iref, 0, 4) ne KEY1 || substr($$iref, 10, 6) ne KEY2);

        #compression block size (bsize-xlen-19)
        $bsize = scalar(unpack('v', substr($$iref, 16, 2))) + 1; #+1 because BSIZE is always too small
        die "ERROR: Invalid BSIZE. Block may be corrupt." if($bsize > 65536);
        last if(length($$iref) < $bsize); #check if string contains full block

        #inflate and concatenate to output
        inflate(\ (substr($$iref, 0, $bsize)), $oref);
    }

    #fail on trailing block if indicated
    die "ERROR: Could not read compression block" if(length($$iref));

    return $oref;
}


#tells the virtual offset of the next alignment record
sub tell_next_alignment {
    my $data = shift; #ref
    my $header = shift;

    return if(!$$data);

    #predeclare values for efficienty
    my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,
	$tlen,$bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar);
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
	$flag = ($flag_nc>>16); #flag is probably not usful here
	$n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	
	#get read name
	$read_name = substr($$data, $substr_off, $l_read_name);
	$substr_off += $l_read_name;
	next unless(substr($read_name, -1) eq "\0"); #last character must be null
	
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
    return EOF_BLOCK
}

sub thread_validate {
    my $infiles = shift;
    my $size = shift;
    my $count = shift;
    my $window = shift;
    my $sections = shift;
    my $O = shift;

    #fix thread signalling
    if(is_thread){
        select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
        select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

        $SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
        $SIG{INT} = sub {exit(130)};
    }

    #open input file
    my $IN = get_handle('<', $file);

    #now read alignment sections
    while(my $section = shift @$sections){
        my $start = ($section-1) * $window;
        my $end   = $section * $window;
        $end = $size if($end > $size); #don't go after end

        #adjust to start of next block
        seek_next_block($IN, $start, 0) if($start > 0);
	my $pos = tell($IN);
	die "ERROR: Cannot find BGZF block within reasonable interval, file may be corrupt\n"
	    if($pos-$start >= 65536);

        #validate all blocks in section
        while(tell($IN) < $end){
            my $block;
	    readblock($IN);
            #inflate(readblock($IN), \$block);
        }

        #validate one more (avoids weirdness with bad blocks at boundaries)
        if(tell($IN) < $size){
            my $block;
	    readblock($IN);
            #inflate(readblock($IN), \$block);
        }
    }

    #clear any open file handles
    close_handles();

    return 1;
}

END {
    unlink(@CLEANUP);
    unlink(@temp_files) unless(is_thread);
}

#the C code to implement MPI calls from perl
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

unsigned long merge_fasthash (unsigned long hash1, unsigned long hash2) {
    uint64_t hash3;
    hash3 = (uint64_t)hash1 + (uint64_t)hash2;

    return hash3;
}

#define bam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const unsigned char seq_nt16_comp_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     8, 4, 2, 1, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 8, 7, 4, 11,15,15, 2, 13,15,15, 3, 15,12,15,15,
    15,15,10, 6,  1,15,14, 9, 15, 5,15,15, 15,15,15,15,
    15, 8, 7, 4, 11,15,15, 2, 13,15,15, 3, 15,12,15,15,
    15,15,10, 6,  1,15,14, 9, 15, 5,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

int w = 0; //global warning flag (print message just once)

int seq_bam2string(SV* ref, int len, int reverse) {
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

int qual_0_to_33(SV* ref, int len, int q64, int reverse) {
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);
    
    //shift quality
    int i;
    if(q64 == 0){ //Illumina+33 as Phred+0 to Phred+33
        for (i = 0; i < len; ++i){
	    qual[i] += 33;
	    if(qual[i] > 75 && w++ < 1) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
        }
    }
    else if(q64 == 1){ //Illumina+64 as Phred+0 to Phred+33
        for (i = 0; i < len; ++i){
	    if(qual[i] < 31) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] += 2;
	    if(qual[i] > 75 && w++ < 1) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
        }
    }
    else if(q64 == 2){ //Solexa+64 as Phred+0 to Phred+33
        int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 26) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-31;
	    qual[i] = 33 + qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
	    if(qual[i] > 75 && w++ < 1) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	}
    }
    else{
        croak("ERROR: Invalid quality encoding conversion value\n");
    }

    if (reverse != 0) {
        uint8_t buf;
        uint8_t *p1 = qual; //start character
	uint8_t *p2 = qual + len - 1; //end character
	while (p1 < p2) {
	    buf = *p1;
	    *p1++ = *p2;
	    *p2-- = buf;
        }
    }
    
    sv_setpvn(scalar, qual, len);
}

//fixes string formated in phred+33
int qual_33_fix(SV* ref, int len, int q64, int reverse) {
    if(q64 == 0 && reverse == 0) return;
    
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);
    
    //shift quality
    int i;
    if(q64 == 0){ //Illumina+33 to Phred+33
        //do nothing
    }
    else if(q64 == 1){ //Illumina+64 to Phred+33
        for (i = 0; i < len; ++i){
	    if(qual[i] < 64) croak("ERROR: Qaulity values are too low to be in Phred+64 format\n");
	    qual[i] -= 31;
	    if(qual[i] > 75 && w++ < 1) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
        }
    }
    else if(q64 == 2){ //Solexa+64 to Phred+33
        int qsol;
        for (i = 0; i < len; ++i){
	    if(qual[i] < 59) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-64;
	    qual[i] = 33 + qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
	    if(qual[i] > 75 && w++ < 1) warn("WARNING: Qaulity values appear too high to be in Phred+33 format\n");
	}
    }
    else{
	croak("ERROR: Invalid quality encoding conversion value\n");
    }
    
    if (reverse != 0) {
	uint8_t buf;
	uint8_t *p1 = qual; //start character
	    uint8_t *p2 = qual + len - 2; //end character
	    while (p1 < p2) {
		buf = *p1;
		*p1++ = *p2;
		*p2-- = buf;
	}
    }

    sv_setpvn(scalar, qual, len);
}

//stats of quality range
int qual_stats(SV* ref, int len, SV* stat_ref) {
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);
    
    HV* stat;
    if (SvTYPE(SvRV(stat_ref))==SVt_PVHV)
        stat = (HV*)SvRV(stat_ref);
    else
        croak("ERROR: Variable is not a hash refrence");
    
    int count = 1;
    int min = qual[0];
    int max = qual[0];
    int sum = qual[0];
    int mean;
    
    int i;
    for (i = 1; i < len; ++i) {
	if(qual[i] < min) min = qual[i];
	if(qual[i] > max) max = qual[i];
	sum += qual[i];
    }
    
    //check stored values
    if(hv_exists(stat, "min", 3)){
	int amin = (int)SvIV(*hv_fetch(stat, "min", 3, 0));
	int amax = (int)SvIV(*hv_fetch(stat, "max", 3, 0));
	int asum = (int)SvIV(*hv_fetch(stat, "sum", 3, 0));
	int alen = (int)SvIV(*hv_fetch(stat, "len", 3, 0));
	int acount = (int)SvIV(*hv_fetch(stat, "count", 5, 0));
	if(amin < min) min = amin;
	if(amax > max) max = amax;
	sum += asum;
	len += alen;
	count += acount;
    }

    //set values
    mean = sum/len;
    hv_store(stat, "min", 3, newSViv(min), 0);
    hv_store(stat, "max", 3, newSViv(max), 0);
    hv_store(stat, "sum", 3, newSViv(sum), 0);
    hv_store(stat, "len", 3, newSViv(len), 0);
    hv_store(stat, "mean", 4, newSViv(mean), 0);
    hv_store(stat, "count", 5, newSViv(count), 0);
}

short reg2bin(int beg, int end){
    if (beg != end) --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
    return 0;
}
