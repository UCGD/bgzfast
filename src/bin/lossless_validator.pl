#!/usr/bin/perl
use forks;
use forks::shared;

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil :sys_wait_h);
use Fcntl;
use Compress::Zlib;
#use Digest::MD5 qw(md5);
BEGIN {
    binmode(STDIN);
    binmode(STDOUT);
    select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
    select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

    $SIG{INT} = sub {exit(130)}; #allow interupt
    mkdir("$ENV{HOME}/.Inline") if(! -d "$ENV{HOME}/.Inline");
}
use Inline qw(C);

#usage statement
my ($exe) = $0 =~ /([^\/]+)$/;
my $usage = "
Usage:

     $exe <bam_file> <fastq_file1> <fastq_file2> ...

     Validates that a BAM file contains the same reads as the given FASTQ files.
     This ensures that the BAM is in fact lossless. The mechanism used is a fast
     checksum ID generator. The checksum for each read is summed, so file order
     is irrelevent.

Options:
     cpus|c     <INT>     CPUs to use for conversion

     restore    <TAG>     Restore original quality string using values from the
                          given BAM tag

     q64                  Fix Phred+64 qaulity scores in fastq (should be Phred+33)

     solexa               Fix Solexa odds based scores in fastq (should be Phred+33)

     help|?               Prints this usage statement

";

#predeclare default option
my $cpus = 1;
my $q64;
my $solexa;
my $fix_qual;
my $restore;

#get options from command line
my @argv = @ARGV; #backup
GetOptions("cpus|c=i" => \$cpus,
           "q64" => \$q64,
           "solexa" => \$solexa,
           "restore=s" => \$restore,
           "help|?" => sub{print $usage; exit(0)});
$fix_qual = 1 if($q64);
$fix_qual = 2 if($solexa);
die "ERROR: The options -q64 and -solexa are mututally exclusive. You cannot use both.\n"
    if($q64 && $solexa);

#get input files
my $bam = shift;
my @fastqs = @ARGV;

#give usage statement
if(!$bam || ! @fastqs){
    print STDERR "WARNING: Failure to supply both BAM and FASTQ for validataion\n";
    print $usage;
    exit(1);
}

#error check options
my $err;
foreach my $file (grep {defined($_)} ($bam, @fastqs)){
    $err .= "ERROR: File does not exist: $file\n" if(! -f $file);
}
die $err if($err);

#make sections of right size
my ($size, $count, $window)= get_chunk_stats($bam);
my @sections = (1..$count);
share(@sections); #share is non destructive for forks::shared

#launch jobs in parallel
my @threads;
my %param = (FIX_QUAL => $fix_qual, RESTORE => $restore);
for(my $i = 1; $i < $cpus; $i++){
    my $thr = threads->new({context => 'list'}, \&thread_run, $bam, \@fastqs, \@sections, \%param);
    push(@threads, $thr);
}
my ($md5a, $md5b) = thread_run($bam, \@fastqs, \@sections, \%param); #let main process join in

#gather results
foreach my $thr (@threads){
    my ($digesta, $digestb) = $thr->join; #wait on threads
    $md5a = merge_fasthash($digesta, $md5a);
    $md5b = merge_fasthash($digestb, $md5b);
}

if($md5a eq $md5b){
    print "SUCCESS: BAM is lossless\n";
}
else{
    print "FAILURE: BAM is not lossless\n";
}

#finished
exit();

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------

#convenience method to see if I am in a thread or not
{my $tid; BEGIN{$tid = threads->tid if exists $INC{'threads.pm'}};
sub is_thread {
    my $is_thread = 1 if(defined($tid) && $tid != threads->tid);
    return $is_thread;
}}

sub get_chunk_stats {
    my $file = shift;

    my $size += (stat($file))[7];
    my $count = ceil($size/33554432); #32Mb sections
    my $window = ceil($size/$count);

    return ($size, $count, $window);
}

#main bam conversion done here by each process
sub thread_run {
    my $bam      = shift;
    my $fastqs   = shift;
    my $sections = shift;
    my $param    = shift;

    #get needed parameters
    my $fix_qual = $param->{FIX_QUAL};
    my $restore  = $param->{RESTORE};
    my ($bsize,  $bcount, $bwindow) = get_chunk_stats($bam);

    #data for fastqs as a group
    my $ftotal = 0;
    my %fstats;
    foreach my $f (@$fastqs){
	my $fsize = (stat($f))[7];
	my $foff_beg = $ftotal;
	$ftotal += $fsize;
	my $foff_end = $foff_beg + $fsize;
	$fstats{$f} = [$fsize, $foff_beg, $foff_end];
    }
    my $fwindow = ceil($ftotal/$bcount);

    #declare variables
    my $header;
    my $header_off;
    my $header_block;

    #open input file
    sysopen(my $BAM, $bam, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $bam: $!";

    #get header
    inflate(readblock($BAM), \$header_block);

    my $magic = substr($header_block, 0, 4);
    $header_off += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
        if($magic ne "BAM\1");

    #grow header block to needed size
    my $l_text = unpack('l<', substr($header_block, $header_off, 4));
    $header_off += 4;
    while((length($header_block)-$header_off < $l_text+4)){
        inflate(readblock($BAM), \$header_block);
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
            inflate(readblock($BAM), \$header_block);
        }
        my $l_name = unpack('l<', substr($header_block, $header_off, 4));
        $header_off += 4;

        while(length($header_block)-$header_off < $l_name+4){ #grow if needed
            inflate(readblock($BAM), \$header_block);
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
    $header_off = systell($BAM); #make offset be position of first block after header

    #now read alignment sections
    my $md5a = 0;
    my $md5b = 0;
    while(my $section = shift @$sections){
        my $start = ($section-1) * $bwindow;
        my $end   = $section * $bwindow;
        $start = $header_off if($start < $header_off);
        $end = $bsize if($end > $bsize); #don't go after end

        #adjust to start of next block
        my $pos = seek_next_block($BAM, $start, 0);

        #read in entire section
	my $data = '';
        my $dlen = ($end-$start) + 2*65536; #data length with overhang of two blocks
	$dlen = ($bsize-$pos) if($pos+$dlen > $bsize); #adjust to end of file
	my $stat = sysread($BAM, $data, $dlen, 0);
	die "ERROR: Failure to read all expected data from file: $bam\n" if($dlen != $stat);

	#make virtual positions
	my $vstart = 0;
	my $vend = $end-$pos;
	my $vpos = $dlen-length($data);

        #get first block in section
        my $block = ($section == 1) ? $header_block : ''; #piece of alignment stuck at end of header
        inflate(substrblock(\$data, \$vpos), \$block);

        #find first alignment
        my $offset = tell_next_alignment(\$block, $header);
        while(! defined($offset) && $vpos < $vstart+65536 && $vpos < $dlen){ #get past short and empty blocks
            inflate(substrblock(\$data, \$vpos), \$block);
            $offset = tell_next_alignment(\$block, $header);
        }
        die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($block));
        substr($block, 0, $offset) = ''; #chop off leading data

        #get all blocks in section
	while($vpos < $vend){
	    inflate(substrblock(\$data, \$vpos), \$block);
	    $md5a = bam2checksum(\$block, $md5a, 0, $restore);
	}

        #get last block and adjust split terminating line
        my $tail = '';
        inflate(substrblock(\$data, \$vpos), \$tail);
        $offset = tell_next_alignment(\$tail, $header);
        while(! defined($offset) && $vpos < $vend+65536 && $vpos < $dlen){ #get past short and empty blocks
            inflate(substrblock(\$data, \$vpos), \$tail);
            $offset = tell_next_alignment(\$tail, $header);
        }
        die "ERROR: No alignments found in section: $section\n" if(! defined($offset) && length($tail));
        $block .= substr($tail, 0, $offset) if($offset); #chop off trailing data

        #get md5 of reads
        $md5a = bam2checksum(\$block, $md5a, 0, $restore);

	#process same section in fastqs
        $start = ($section-1) * $fwindow;
        $end   = $section * $fwindow;
        $end = $ftotal if($end > $ftotal); #don't go after end

	#files to use for this section
	my @files = grep {overlap($start, $end, $fstats{$_}[1], $fstats{$_}[2])} @$fastqs;
	foreach my $f (@files){
	    #get offsets in combined file space
	    my ($fsize, $foff_beg, $foff_end) = @{$fstats{$f}};
	    my $adjust = $foff_beg;
	    $foff_beg = $start if($start > $foff_beg);
	    $foff_end = $end if($end < $foff_end);

	    #open file and adjust to start and end of fastq entries
	    sysopen(my $FQ, $f, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $f: $!";
	    $foff_end = seek_next_fq($FQ, $foff_end-$adjust);
	    $foff_beg = seek_next_fq($FQ, $foff_beg-$adjust);

	    #make md5 for everything up until end of file or section
	    my $data = '';
	    my $dlen = $foff_end-$foff_beg;
	    my $stat = sysread($FQ, $data, $dlen, 0); #max block 
	    die "ERROR: Failure to read all expected data from file: $f\n" if($dlen != $stat);

	    my $vpos = 0;
	    while(my $read = substrfq(\$data, \$vpos,$fix_qual)){
		my $digest = fasthash($read);
		$md5b = merge_fasthash($digest, $md5b);
	    }
	}
    }

    return ($md5a, $md5b);
}

sub overlap {
    my $Abeg = shift;
    my $Aend = shift;
    my $Bbeg = shift;
    my $Bend = shift;

    return 0 if($Aend <= $Bbeg);
    return 0 if($Abeg >= $Bend);
    return 1;
}

{#predeclare values for efficienty
my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,$tlen,
    $bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,$seq,@qual,$qual);
my ($tag, $val_type, $value, $val_len, $sub_type, $sub_len);
my ($cl_ref,$cl_query,$match);
my ($substr_off, $data_len, $block_offset);
my ($data, $md5, $fix_qual, $restore);

my @seq_index;
BEGIN {@seq_index = qw(= A C M G R S V T W Y H K D B N);} #force it to exist

sub bam2checksum {
    $data = shift; #ref
    $md5  = shift;
    $fix_qual  = shift || 0;
    $restore  = shift;

    #collect and process bam alignments
    $data_len = length($$data);
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
        $read_name = substr($$data, $substr_off, $l_read_name);
        $substr_off += $l_read_name;
        chop($read_name); #removes trailing null

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

        #process seq string (slow)
        my $count = int(($l_seq+1)/2);
        $seq = substr($$data, $substr_off, $count);
        $substr_off += $count;
        convert_seq(\$seq, $l_seq, ($flag & 16));

        #process quality string (slow)
        my $qual = substr($$data, $substr_off, $l_seq);
        $substr_off += $l_seq;
	if(!$restore){
	    convert_qual(\$qual, $l_seq, ($flag & 16), $fix_qual);
	}
        else{ #process tags to get original quality values
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

                if($tag eq $restore){
                    die "ERROR: Wrong datatype for restored qaulity value\n" if($val_type ne 'Z');
                    die "ERROR: Restored quality value does not match sequence length\n" if($val_len-1 != $l_seq);
                    $value = substr($$data, $substr_off, $val_len-1); #-1 to ignore null
                    $substr_off += $val_len;
                    last; #short circuit loop
                }
                else{ #skip past value since it is not the right one
                    $substr_off += $val_len;
                }
            }

            #restore value
            if(!$value){
                warn "WARNING: Original quality value not found for $read_name\n";
		convert_qual(\$qual, $l_seq, ($flag & 16), $fix_qual);
            }
            else{
                fix_qual(\$value, $l_seq, ($flag & 16), $fix_qual) if(($flag & 16) || $fix_qual);
                $qual = $value;
            }
        }

	#get checksum of single read entry and merge it with existing entries (ignore pair in ID)
	my $digest = fasthash("\@$read_name\n$seq\n+\n$qual\n");
	$md5 = merge_fasthash($digest, $md5);
    }

    $$data = substr($$data, $block_offset); #remove what I processed off of the block

    return $md5;
}}

#not a perfect merge but differences are sytematic
sub merge_md5 {
    my $md5a = shift;
    my $md5b = shift;

    return $md5a if(!defined($md5b));
    return $md5b if(!defined($md5a));

    my @A = unpack('LLLL', $md5a);
    my @B = unpack('LLLL', $md5b);
    return pack('LLLL', $A[0]+$B[0], $A[1]+$B[1], $A[2]+$B[2], $A[3]+$B[3]);
}

sub systell {
    return sysseek($_[0], 0, 1);
}

sub seek_next_fq {
    my $FH = shift; #filehandle
    my $start = shift; #where to start looking from
    my $wence = shift || 0;

    #go to position
    seek($FH, $start, $wence); #regular seek since readline is buffered IO
    my $pos = tell($FH); #regular tell since readline is buffered IO
    my $size = (stat($FH))[7]; #file size

    return sysseek($FH, 0, 0) if($pos == 0);
    return sysseek($FH, $size, 0) if($pos == $size);

    #get 4 line
    $pos += length(readline($FH)); #skip partial lines
    my @l = grep {defined($_)} map {scalar(readline($FH))} (1..4); #4 lines
    my @c = map {substr($_, 0, 1)} @l; #1st char of each line
    if(@l != 4){
	return sysseek($FH, $size, 0) if(tell($FH) == $size); #end of file
	die "ERROR: Failure to read sufficient data from file\n";
    }

    #use pattern to adjust
    if(   $c[0] =~ /^\@/ && $c[1] =~ /^[ATCGURYKMSWBDHVN]/ &&
	  $c[2] =~ /^\+/ && 33 <= ord($c[3]) && ord($c[3]) <= 126){
	#do nothing
    }
    elsif($c[1] =~ /^\@/ && $c[2] =~ /^[ATCGURYKMSWBDHVN]/ &&
	  $c[3] =~ /^\+/ && 33 <= ord($c[0]) && ord($c[0]) <= 126){
	$pos += length(shift @l);
    }
    elsif($c[2] =~ /^\@/ && $c[3] =~ /^[ATCGURYKMSWBDHVN]/ &&
	  $c[0] =~ /^\+/ && 33 <= ord($c[1]) && ord($c[1]) <= 126){
	
	$pos += length(shift @l);
	$pos += length(shift @l);
    }
    elsif($c[3] =~ /^\@/ && $c[0] =~ /^[ATCGURYKMSWBDHVN]/ &&
	  $c[1] =~ /^\+/ && 33 <= ord($c[2]) && ord($c[2]) <= 126){
	
	$pos += length(shift @l);
	$pos += length(shift @l);
	$pos += length(shift @l);
    }
    else{
	die "ERROR: Pattern does not match expected format for FASTQ\n"
    }

    return sysseek($FH, $pos, 0);
}

sub seek_next_block {
    my $FH = shift; #filehandle
    my $start = shift; #where to start looking from
    my $wence = shift || 0;

    #go to position
    sysseek($FH, $start, $wence);
    my $pos = systell($FH);
    my $size = (stat($FH))[7]; #file size

    #not enough space left for a block
    return sysseek($FH, $size, 0) if($size-$start < 28);

    #read in some data to check
    my $data;
    sysread($FH, $data, 65536, 0); #max block

    #find block keys in data
    my $offset = 0;
    my $key1 = pack('H8', '1f8b0804'); #4 byte static header key
    my $key2 = pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes)
    while(1){
        #match first static key
        $offset = index($data, $key1, $offset);
        return sysseek($FH, $size, 0) if($offset == -1);

        #match second static key
        my $offset2 = index($data, $key2, $offset);
        return sysseek($FH, $size, 0) if($offset2 == -1);

        #second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
        if($offset2-$offset == 10){
            last;
        }
        else{
            $offset += 4; #jump over key1 match and try again
        }
    }
    $pos += $offset; #adjust file position with offset

    return sysseek($FH, $pos, 0);
}

#assumes file position is set to start of block
sub readblock {
    my $FH = shift;

    #read BGZF header
    my $data;
    my $stat = sysread($FH, $data, 18, 0);
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
    $stat = sysread($FH, $data, $c_data_size + 8, length($data)); #the +8 is for the footer
    die "ERROR: Could not read compression block\n"
        if($stat != $c_data_size + 8);

    return \$data;
}

sub substrfq {
    my $str_ref  = shift;
    my $pos_ref  = shift;
    my $fix_qual = shift;

    return undef if($$pos_ref >= length($$str_ref));

    my $read;
    for(my $i = 0; $i < 4; $i++){
	my $j = index($$str_ref, "\n", $$pos_ref);
	if($j == -1){
	    $read .= substr($$str_ref, $$pos_ref);
	    $$pos_ref = length($$str_ref);
	    die "ERROR: Malformed FASTQ format" if($i != 3);
	}

	if($i == 0){
	    $read = substr($$str_ref, $$pos_ref, $j+1-$$pos_ref);
	    $$pos_ref = $j+1;
	    $read =~ s/[\s\t]+[^\n]+$//; #remove comments
	    $read =~ s/\/\d+$//; #remove read pairing tags
	}
	elsif($i == 1){
	    $read .= substr($$str_ref, $$pos_ref, $j+1-$$pos_ref);
	    $$pos_ref = $j+1;
	}
	elsif($i == 2){
	    $read .= "+\n";
	    $$pos_ref = $j+1;
	}
	elsif($i == 3){
	    if($fix_qual){ #fix FASTQ quality values
		my $qual = substr($$str_ref, $$pos_ref, $j+1-$$pos_ref);
		$$pos_ref = $j+1;
		fix_qual(\$qual, length($qual), 0, $fix_qual);
		$read .= $qual;
	    }
	    else{
		$read .= substr($$str_ref, $$pos_ref, $j+1-$$pos_ref);
		$$pos_ref = $j+1;
	    }
	}
    }

    return $read;
}

sub substrblock {
    my $str_ref = shift;
    my $pos_ref = shift;

    #read BGZF header
    my $data = substr($$str_ref, $$pos_ref, 18);
    my $stat = length($data);
    $$pos_ref += $stat;
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
    $data .= substr($$str_ref, $$pos_ref, $c_data_size + 8); #the +8 is for the footer
    $stat = length($data) - 18;
    $$pos_ref += $stat;
    die "ERROR: Could not read compression block\n"
        if($stat != $c_data_size + 8);

    return \$data;
}

sub inflate {
    my $buffer = shift; #reference
    my $ref    = shift; #reference

    return if(!$buffer);

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 131072);
    die "ERROR: Failed to create zlib inflation object with status: $stat\n" if($stat);

    substr($$buffer, 0, 18, '');
    $$ref .= scalar($i_obj->inflate($buffer));

    #validate the length and crc32
    die "ERROR: Trailing garbage in compression block\n" if(length($$buffer) != 8);
    
    my ($crc32, $isize) = unpack('L<L<', $$buffer);
    if($isize != $i_obj->total_out()){ #size does not match
	die "ERROR: The expected ISIZE of the uncompressed block does not match\n";
    }
    if($crc32 != crc32(substr($$ref, -$isize, $isize))){ #crc32 does not match
	die "ERROR: The expected CRC32 of the uncompressed block does not match\n";
    }

    return;
}

sub deflate {
    my $buffer = shift; #reference
    my $ref = shift; #reference

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 131072, -Level => 4);
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

    return if(!$$data);

    #predeclare values for efficienty
    my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,
        $tlen,$bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,@qual);
    my ($substr_off, $l_ref, $match, $cl_ref, $cl_query);

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

        return $i; #seems to be a valid bam line
    }

    return;
}

sub eof_block {
    return pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000');
}

#------------------------------------------------------------------------------
#------------------------------- Inline C Code --------------------------------
#------------------------------------------------------------------------------

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
    
    uint64_t hash;
    hash = 0;

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

int fix_qual(SV* ref, int len, int reverse, int fix_flag) {
    //note these qualities are already expected to be phred+33
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    //char* buf = (char*)malloc((len+1)*sizeof(char));
    //buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(fix_flag == 2){ //recalculate solexa encoding
        int qsol;
	for (i = 0; i < len; ++i){
	    if(qual[i] < 59) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
	    qsol = qual[i]-64;
	    qual[i] = 33 + qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
	}
    }
    else if(fix_flag != 0){ //other illumina encoding
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

int convert_qual(SV* ref, int len, int reverse, int fix_flag) {
    SV* scalar = SvRV(ref);
    uint8_t* qual = (uint8_t*)SvPV_nolen(scalar);

    char* buf = (char*)malloc((len+1)*sizeof(char));
    buf[len] = 0;

    //recalculate quality shift of 64
    int i;
    if(fix_flag == 2){ //recalculate solexa encoding
        int qsol;
        for (i = 0; i < len; ++i){
            if(qual[i] < 26) croak("ERROR: Qaulity values are too low to be in Solexa format\n");
            qsol = qual[i]-31;
            qual[i] = qsol + 10 * log10(1 + pow(10, (double)qsol/-10)); //change scaling equation
        }
    }
    else if(fix_flag != 0){ //other illumina encoding
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
