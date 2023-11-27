#------------------------------------------------------------------------
#----                   BGZFast::bam_manipulator                     ----
#------------------------------------------------------------------------
package BGZFast::bam_manipulator;

use strict;
use warnings;
use Compress::Zlib;
use Fcntl qw(:seek O_RDONLY); #SEEK_SET, SEEK_CUR, SEEK_END

use base qw(BGZFast::seqorder_manipulator);

our $VERSION = '0.01';

#------------------------------------------------------------------------
#--------------       OBJECT INITIALIZATION METHODS        --------------
#------------------------------------------------------------------------

sub new {
    return shift->SUPER::new(@_);
}

sub init {
    return shift->SUPER::init(@_);
}

sub bai {return shift->idx;}
sub idx {
    my $self = shift;

    return $self->{BAI} if($self->{BAI});

    my $file = $self->__file->{PATH};
    (my $base = $file) =~ s/\.bam$//;
    my $bai_file = (!-f "$file.bai" && -f "$base.bai") ? "$base.bai" : "$file.bai";
    if(! -f $bai_file){
	if(my $exe = File::Which::which('sambamba') || File::Which::which('samtools')){
	    my $stat = system($exe => 'index' => $file);
	    unlink($bai_file) if($stat); #on failure
	}
	die "ERROR: ".__PACKAGE__." could not generate BAM index files\n" if(!-f $bai_file);
    }

    $self->__parse_bai($bai_file);

    return $self->{BAI};
}

sub next_alignment {
    my $self = shift;
    my $buffer = $self->__buffer;
    my $align = $self->__alignment;

    #process bam alignment
    my $data_len = $buffer->{LENGTH};
    my $block_offset = $buffer->{POS};
    while(1){
	while($block_offset >= $data_len - 4){ #load new block if near end
	    return unless($self->__read_block);
	    my $block = $self->__block;
	    my $e_len = $buffer->{LENGTH} + $block->{I_LENGTH}; #expected length
	    $self->__inflate_to_buffer;
	    $data_len = $buffer->{LENGTH};
	    $block_offset -= $e_len - $data_len; #accounts for potential truncation
	}
        my $substr_off = $block_offset;

        my ($block_size) = unpack('l<', substr($buffer->{DATA}, $substr_off, 4));
        $substr_off += 4;
        while($substr_off+$block_size > $data_len){ #current alignment is spit across next block
	    if(!$self->__read_block){
		$self->{ALIGNMENT} = undef;
		return;
	    }

	    my $block = $self->__block;
	    my $e_len = $buffer->{LENGTH} + $block->{I_LENGTH}; #expected length
	    $self->__inflate_to_buffer;
	    $data_len = $buffer->{LENGTH};
	    $block_offset -= $e_len - $data_len; #accounts for potential truncation
	    $substr_off = $block_offset + 4;
	}

	#get entire alignment
	$align->{DATA} = substr($buffer->{DATA}, $substr_off-4, $block_size+4);
        $block_offset = $substr_off+$block_size; #end of current alignment

	#get position part of alignment entry (ignore the rest)
        ($align->{REFID},$align->{POS},undef,undef,undef,$align->{NEXT_REFID},$align->{NEXT_POS})
            = unpack('l<l<L<L<l<l<l<', substr($align->{DATA}, 4, 28));

	#set file/buffer position info
	($align->{FILE_OFFSET}, $align->{BLOCK_OFFSET}) = $self->block_tell;
	$align->{LENGTH} = $block_size+4; #length of total alignment entry
	$buffer->{POS} += $align->{LENGTH};

	$self->{ALIGNMENT} = $align;
	return $align;
    }    
}

#parse header from bam
sub __parse_header {
    my $self = shift;
    my $IN = $self->__handle;	
    my $header = $self->__header;

    #backup current position then move to start of file
    my $c_file_pos = $self->__file->{POS};
    my $c_real_pos = sysseek($IN, 0, SEEK_CUR);
    sysseek($IN, 0, SEEK_SET);

    #read in the header block
    my $h_buffer = '';
    my $h_block = {};
    $self->__read_block($h_block); #reads into $h_block
    $header->{D_LENGTH} += $h_block->{D_LENGTH};

    #inflate header
    $self->__inflate(\ ($h_block->{DATA}), \$h_buffer);

    my $h_offset = 0;
    my $magic = substr($h_buffer, $h_offset, 4);
    $h_offset += 4;
    die "ERROR: Magic string mismatch. This does not appear to be a bam file\n"
        if($magic ne "BAM\1");

    #grow header block to needed size
    my $l_text = unpack('l<', substr($h_buffer, $h_offset, 4));
    $h_offset += 4;
    while((length($h_buffer)-$h_offset < $l_text+4)){
	$self->__read_block($h_block);
	$header->{D_LENGTH} += $h_block->{D_LENGTH};
	$self->__inflate(\ ($h_block->{DATA}), \$h_buffer);
    }

    #get header text
    my $text = substr($h_buffer, $h_offset, $l_text);
    $h_offset += $l_text;
    $text =~ s/\0$//g; #remove null padding
    $header->{TEXT} = $text;

    #get header reference count
    my $n_ref = unpack('l<', substr($h_buffer, $h_offset, 4));
    $h_offset += 4;
    $header->{N_REF} = $n_ref;

    #get each reference
    for(my $i = 0; $i < $n_ref; $i++){
        while(length($h_buffer)-$h_offset < 4){ #grow if needed
	    $self->__read_block($h_block);
	    $header->{D_LENGTH} += $h_block->{D_LENGTH};
	    $self->__inflate(\ ($h_block->{DATA}), \$h_buffer);
        }
        my $l_name = unpack('l<', substr($h_buffer, $h_offset, 4));
        $h_offset += 4;

        while(length($h_buffer)-$h_offset < $l_name+4){ #grow if needed
	    $self->__read_block($h_block);
	    $header->{D_LENGTH} += $h_block->{D_LENGTH};
	    $self->__inflate(\ ($h_block->{DATA}), \$h_buffer);
        }
        my $name = substr($h_buffer, $h_offset, $l_name);
        $h_offset += $l_name;
        $name =~ s/\0$//g; #remove null padding
        $header->{REF}[$i]{NAME} = $name;

        my $l_ref = unpack('l<', substr($h_buffer, $h_offset, 4));
        $h_offset += 4;
        $header->{REF}[$i]{LENGTH} = $l_ref;

	#add conversion
	$header->{NAME2SEQID}{$name} = $i;
	$header->{NAME2LENGTH}{$name} = $l_ref;
    }
    $header->{NAME2SEQID}{'*'} = -1; #unmapped chr
    $header->{NAME2LENGTH}{'*'} = 0; #unmapped chr

    #separate header from alignment tail
    my $h_data = substr($h_buffer,0,$h_offset,'');

    #identify real and virtual offset for end of header (start of alignments)
    $header->{TAIL_FILE_OFFSET} = $h_block->{FILE_OFFSET}; #real offset
    $header->{TAIL_BLOCK_OFFSET} = $h_block->{I_LENGTH} - length($h_buffer); #virtual offset

    #deflate header (and tail)
    $self->__deflate(\ (substr($h_data, 0, 65280, '')), \ ($header->{DATA})) while(length($h_data));
    $self->__deflate(\ (substr($h_buffer, 0, 65280, '')), \ ($header->{TAIL})) while(length($h_buffer));

    #restore file positions
    $self->__file->{POS} = $c_file_pos;
    sysseek($IN, $c_real_pos, SEEK_SET);

    return $header;
}

sub __parse_bai {
    my $self = shift;
    my $file = shift;

    #read in entire file
    my @bai; #list representation
    my $data; #raw data
    my $size = (stat($file))[7];
    sysopen(my $IN, $file, O_RDONLY) or die "ERROR: Could not open file $file: $!";
    binmode($IN);
    sysread($IN, $data, $size);
    close($IN);

    #check magic string
    my $str_off = 0; #offset in substring
    if(substr($data, $str_off, 4) ne "BAI\1"){
	die "ERROR: File does not appear to be a BAM index file\n"
    }
    $str_off += 4;
    
    #get number of reference sequences
    my $n_ref = unpack('l<', substr($data, $str_off, 4));
    $str_off += 4;

    #list of indices (n=n_ref)
    for(my $i = 0; $i < $n_ref; $i++){
	$bai[$i] = {BINS => [], LINEAR => [], PSEUDO => undef}; #initialize index of each contig

	#get number of distinct bins (for the binning index)
	my $n_bin = unpack('l<', substr($data, $str_off, 4));
	$str_off += 4;
	
	#list of distinct bins (n=n_bin)
	my $bai_b = $bai[$i]{BINS};
	for(my $j = 0; $j < $n_bin; $j++){
	    #get distinct bin & number of chunks
	    my ($bin, $n_chunk) = unpack('L<l<', substr($data, $str_off, 8));
	    $str_off += 8;

	    #fill in pseudo-bin
	    if($bin == 37450){
		my $bai_p = $bai[$i]{PSEUDO} = []; #initialize chunk list for pseudo bin
		die "ERROR: Corrupt BAI index for pseudo bin\n" if($n_chunk != 2); #always 2

		#get (virtual) file offset of the start and end of placed unmapped reads
		my ($unmapped_beg, $unmapped_end) = unpack('Q<Q<', substr($data, $str_off, 16));
		$str_off += 16;

		#offsets are encoded as --> coffset<<16|uoffset
		my $coffset_beg = ($unmapped_beg >> 16);
		my $uoffset_beg = $unmapped_beg & 65535; #mask is (1<<16)-1
		my $coffset_end = ($unmapped_end >> 16);
		my $uoffset_end = $unmapped_end & 65535; #mask is (1<<16)-1

		push(@{$bai_p}, [$coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end]);

		#get number of mapped/unmapped read-segments for this reference
		my ($n_mapped, $n_unmapped) = unpack('Q<Q<', substr($data, $str_off, 16));
		$str_off += 16;

		push(@{$bai_p}, [$n_mapped, $n_unmapped]);

		next;
	    }

	    #list of chunks (n=n_chunk)
	    $bai_b->[$bin] = []; #initialize chunk list for each bin
	    for(my $k = 0; $k < $n_chunk; $k++){
		#get (virtual) file offset of the start and end of the chunk
		my ($chunk_beg, $chunk_end) = unpack('Q<Q<', substr($data, $str_off, 16));
		$str_off += 16;
		
		#chunk offsets are encoded as --> coffset<<16|uoffset
		my $coffset_beg = ($chunk_beg >> 16);
		my $uoffset_beg = $chunk_beg & 65535; #mask is (1<<16)-1
		my $coffset_end = ($chunk_end >> 16);
		my $uoffset_end = $chunk_end & 65535; #mask is (1<<16)-1

		push(@{$bai_b->[$bin]}, [$coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end]);
	    }
	}
	
	#get number of 16kbp intervals (for the linear index)
	my $n_intv = unpack('l<', substr($data, $str_off, 4));
	$str_off += 4;
	
	#list of intervals (n=n_intv)
	my $bai_l = $bai[$i]{LINEAR};
	for(my $j = 0; $j < $n_intv; $j++){
	    #get (virtual) file offset of the first alignment in the interval
	    my $ioffset = unpack('Q<', substr($data, $str_off, 8));
	    $str_off += 8;
	    
	    #interval offsets are encoded as --> coffset<<16|uoffset
	    my $i_coffset = ($ioffset >> 16);
	    my $i_uoffset = $ioffset & 65535; #mask is (1<<16)-1

	    $bai_l->[$j] = [$i_coffset, $i_uoffset];
	}
    }
    
    #get number of unplaced unmapped reads (RNAME *)
    my $n_no_coor = unpack('Q<', substr($data, $str_off, 8));
    $str_off += 8;

    #hash representation
    my $bai = $self->__bai;
    $bai->{PATH} = $file;
    $bai->{ABS_PATH} = Cwd::abs_path($file);
    $bai->{LENGTH} = $size;
    $bai->{INDICES} = \@bai;
    $bai->{N_NO_COOR} = $n_no_coor;

    return $n_ref;
}

sub __bai2string {
    my $self = shift;
    my $bai = $self->bai;

    #add magic
    my $data = "BAI\1";

    #number of reference sequences
    my $n_ref = scalar(@{$bai->{REF}});
    $data .= pack('l<', $n_ref);

    #list of indices (n=n_ref)
    for(my $i = 0; $i < $n_ref; $i++){
	#number of distinct bins (for the binning index)
	my $bai_b = $bai->{INDICES}[$i]{BINS}; #real bins
	my $bai_p = $bai->{INDICES}[$i]{PSEUDO}; #pseudo bin
	my @bins = grep {$bai_b->[$_]} (0..$#$bai_b);
	push(@bins, 37450) if($bai_p); #optional pseudo bin
	my $n_bin = scalar(@bins);
	$data .= pack('l<', $n_bin);
	
	#list of distinct bins (n=n_bin)
	for(my $j = 0; $j < $n_bin; $j++){
	    #distinct bin
	    my $bin = $bins[$j];

	    #fill in pseudo-bin
	    if($bin == 37450){
		#number of chunks
		my $n_chunk = 2; #always 2 for pseudo bin
		$data .= pack('L<l<', $bin, $n_chunk);

		#offsets are encoded as --> coffset<<16|uoffset
		my ($coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end) = @{$bai_p->[0]};
		my $unmapped_beg = ($coffset_beg << 16)|$uoffset_beg;
		my $unmapped_end = ($coffset_end << 16)|$uoffset_end;

		#(virtual) file offset of the start and end of placed unmapped reads
		$data .= pack('Q<Q<', $unmapped_beg, $unmapped_end);

		#number of mapped/unmapped read-segments for this reference
		my ($n_mapped,$n_unmapped) = @{$bai_p->[1]};
		$data .= pack('Q<Q<', $n_mapped, $n_unmapped);

		next;
	    }

	    #list of chunks (n=n_chunk)
	    my $n_chunk = scalar(@{$bai_b->[$bin]});
	    $data .= pack('L<l<', $bin, $n_chunk);
	    for(my $k = 0; $k < $n_chunk; $k++){
		#chunk offsets are encoded as --> coffset<<16|uoffset
		my ($coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end) = @{$bai_b->[$bin][$k]};
		my $chunk_beg = ($coffset_beg << 16)|$uoffset_beg;
		my $chunk_end = ($coffset_end << 16)|$uoffset_end;

		#(virtual) file offset of the start and end of the chunk
		$data .= pack('Q<Q<', $chunk_beg, $chunk_end);
	    }
	}
	
	#number of 16kbp intervals (for the linear index)
	my $bai_l = $bai->{INDICES}[$i]{LINEAR};
	my $n_intv = scalar(@$bai_l);
	$data .= pack('l<', $n_intv);
	
	#list of intervals (n=n_intv)
	for(my $j = 0; $j < $n_intv; $j++){
	    #interval offsets are encoded as --> coffset<<16|uoffset
	    my ($i_coffset, $i_uoffset) = @{$bai_l->[$j]};
	    my $ioffset = ($i_coffset << 16)|$i_uoffset;

	    #(virtual) file offset of the first alignment in the interval
	    $data .= pack('Q<', $ioffset);
	}
    }
    
    #number of unplaced unmapped reads (RNAME *)
    my $n_no_coor = $bai->{N_NO_COOR};
    $data .= pack('Q<', $n_no_coor);

    return $data;
}

sub seqid2name {
    my $self = shift;
    my $id = shift;
    my $header = $self->header;

    if(! exists($header->{REF}[$id])){
	warn "WARNING: No contig exists with seqID $id\n";
	return undef;
    }

    return $header->{REF}[$id]{NAME};
}

sub seqid2seqlength {
    my $self = shift;
    my $id = shift;
    my $header = $self->header;

    if(! exists($header->{REF}[$id])){
	warn "WARNING: No contig exists with seqID $id\n";
	return undef;
    }

    return $header->{REF}[$id]{LENGTH};
}

sub name2seqid {
    my $self = shift;
    my $name = shift;
    my $header = $self->header;

    if(! exists($header->{NAME2SEQID}{$name})){
	warn "WARNING: No contig exists with name $name\n";
	return undef;
    }

    return $header->{NAME2SEQID}{$name};
}

sub name2seqlength {
    my $self = shift;
    my $name = shift;
    my $header = $self->header;

    if(! exists($header->{NAME2LENGTH}{$name})){
	warn "WARNING: No contig exists with name $name\n";
	return undef;
    }

    return $header->{NAME2LENGTH}{$name};
}

#inflate compression block
sub __inflate {
    my $self = shift;
    my $in_buffer  = shift; #reference
    my $out_buffer = shift; #reference

    return if(!$in_buffer);

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 131072);
    die "ERROR: Failed to create zlib inflation object with status: $stat\n" if($stat);

    substr($$in_buffer, 0, 18, '');
    $$out_buffer .= scalar($i_obj->inflate($in_buffer));

    #validate the length and crc32
    die "ERROR: Trailing garbage in compression block\n" if(length($$in_buffer) != 8);
    
    my ($crc32, $isize) = unpack('L<L<', $$in_buffer);
    $$in_buffer = ''; #should now be empty
    
    if($isize != $i_obj->total_out()){ #size does not match
	die "ERROR: The expected ISIZE of the uncompressed block does not match\n";
    }
    if($crc32 != crc32(substr($$out_buffer, -$isize, $isize))){ #crc32 does not match
	die "ERROR: The expected CRC32 of the uncompressed block does not match\n";
    }

    return;
}

#deflate compresssion block
sub __deflate {
    my $self = shift;
    my $in_buffer = shift; #reference
    my $out_buffer = shift; #reference

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 131072, -Level => 4);
    die "ERROR: Failed to create zlib deflation object with status: $stat\n" if($stat);

    my $offset = length($$out_buffer);
    $$out_buffer .= pack('H36', '1f8b08040000000000ff0600424302000000'); #header (BSIZE=0)
    $$out_buffer .= scalar($d_obj->deflate($in_buffer)).scalar($d_obj->flush); #compressed data
    $$out_buffer .= pack('V', crc32($in_buffer)); #CRC32
    $$out_buffer .= pack('V', $d_obj->total_in()); #ISIZE
    substr($$out_buffer, $offset+16, 2, pack('v', $d_obj->total_out+6+19)); #set final BSIZE

    return;
}

#assumes file position is set to start of block
sub __read_block {
    my $self = shift;
    my $block = shift || $self->__block;
    my $IN = $self->__handle;

    #read BGZF header
    my $stat = sysread($IN, $block->{DATA}, 18, 0);
    return undef if($stat == 0); #end of file
    die "ERROR: Failure to read BAM stream\n" if($stat != 18);

    #validate header
    my ($id1,$id2,$cm,$flg,$mtime,$xfl,$os,
        $xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $block->{DATA});
    die "ERROR: Does not appear to be a BGZF file\n"
        if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
           $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);

    #read compression block and footer
    my $c_data_size = $bsize-$xlen-19; #compression block size
    $stat = sysread($IN, $block->{DATA}, $c_data_size + 8, 18); #the +8 is for the footer
    die "ERROR: Could not read compression block\n"
        if($stat != $c_data_size + 8);

    #fill in block values
    my $file = $self->__file;
    $block->{FILE_OFFSET} = $file->{POS};
    $block->{D_LENGTH} = $c_data_size + 18 + 8;
    ($block->{I_LENGTH}) = unpack('L<', substr($block->{DATA}, -4));
    $file->{POS} += $block->{D_LENGTH}; #move file position
    
    return 1;
}

#walks across file looking for the next block header
sub __adjust_to_next_block {
    my $self = shift;
    my $IN = $self->__handle; #filehandle

    #go to position
    my $pos = $self->__file->{POS};
    my $size = $self->__file->{LENGTH}; #file size

    #not enough space left for a block
    return $self->file_seek($size) if($size-$pos < 28);

    #read in some data to check
    my $data;
    sysread($IN, $data, 65536, 0); #max block

    #find block keys in data
    my $offset = 0;
    my $key1 = pack('H8', '1f8b0804'); #4 byte static header key
    my $key2 = pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes)
    while(1){
        #match first static key
        $offset = index($data, $key1, $offset);
        return $self->file_seek($size) if($offset == -1);

        #match second static key
        my $offset2 = index($data, $key2, $offset);
        return $self->file_seek($size) if($offset2 == -1);

        #second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
        if($offset2-$offset == 10){
            last;
        }
        else{
            $offset += 4; #jump over key1 match and try again
        }
    }
    $pos += $offset; #adjust file position with offset

    return $self->file_seek($pos);
}

#walks across uncompressed block looking for alignment header
sub __adjust_to_next_alignment {
    my $self = shift;
    my $buffer = $self->__buffer;
    my $header = $self->header;

    #reset value
    $self->{ALIGNMENT} = undef if($self->{ALIGNMENT});

    #load any deflated blocks into buffer
    $self->__inflate_to_buffer if($self->{BLOCK});

    #get block if buffer is empty
    if(!$buffer->{DATA}){
	return unless($self->__read_block);
	$self->__inflate_to_buffer;
    }

    #keep loading more blocks if buffer is too small
    while($buffer->{LENGTH} - $buffer->{POS} < 65536){
	last unless($self->__read_block);
	$self->__inflate_to_buffer;
    }

    #predeclare values for efficienty
    my ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,$next_refID,$next_pos,
        $tlen,$bin,$mapq,$l_read_name,$flag,$n_cigar_op,$read_name,@cigar,@seq,@qual);
    my ($substr_off, $l_ref, $match, $cl_ref, $cl_query);
    my $n_ref = $header->{N_REF};

    #find correct virtual offset
    L1: for(my $i = $buffer->{POS}; $i < $buffer->{LENGTH}-44; $i++){
        $substr_off = $i; #shift 1 byte to the right each time through loop

        ($block_size,$refID,$pos,$bin_mq_nl,$flag_nc,$l_seq,
         $next_refID,$next_pos,$tlen) = unpack('l<l<l<L<L<l<l<l<l<', substr($buffer->{DATA}, $substr_off, 36));
        $substr_off += 36;

        next unless(44 <= $block_size && $block_size < 65536); #reasonable assumption?
        next unless( 1 <= $l_seq && $l_seq < 10000); #reasonable assumption?

        #bin_mq_nl processing into sub values
        $bin = ($bin_mq_nl>>16);
        $mapq = (($bin_mq_nl>>8) & 255); #mask is (1<<8)-1
        $l_read_name = ($bin_mq_nl & 255); #mask is (1<<8)-1
	next unless($l_read_name < 1000); #reasonable assumption?
        next unless($bin <= 37448); #BAM spec allows for 37448 max bin with 37450 as pseudobin

        #flag_nc processing into sub values
        $flag = ($flag_nc>>16);
	next unless($flag < 4096); #only first 12 bits are used
        $n_cigar_op = ($flag_nc & 65535); #mask is (1<<16)-1
	next unless($n_cigar_op <= $l_seq); #reasonable assumption?

	#grow buffer if necessary
	while($buffer->{LENGTH} < $substr_off+$l_read_name){
	    next L1 unless($self->__read_block);
	    $self->__inflate_to_buffer;
	}

        #get read name
        $read_name = substr($buffer->{DATA}, $substr_off, $l_read_name);
        $substr_off += $l_read_name;
	next unless($read_name =~ /^(?:\*|[!-()+-<>-~][!-~]*)\0$/); #restricted character space

	#grow buffer if necessary
	while($buffer->{LENGTH} < $substr_off+($n_cigar_op*4)){
	    next L1 unless($self->__read_block);
	    $self->__inflate_to_buffer;
	}

        #process cigar into sub values
        @cigar = unpack("L<"x$n_cigar_op, substr($buffer->{DATA}, $substr_off, $n_cigar_op*4));
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

	#ignore TLEN (removed because TLEN is niether consistent nor calculable without both reads)
        #next unless($refID != $next_refID || $pos < 0 || $next_pos < 0 || $next_pos-($pos+$cl_ref) == $tlen);

        #additionally validate contig and positions info (overkill?)
	next unless(-1<= $refID && $refID < $n_ref);
	$l_ref = ($refID != -1) ? $header->{REF}[$refID]{LENGTH} : 0;
	next unless(-1 <= $pos && $pos < $l_ref);
	next unless(-1 <= $next_refID && $next_refID < $n_ref);
	$l_ref = ($next_refID != -1) ? $header->{REF}[$next_refID]{LENGTH} : 0;
	next unless(-1 <= $next_pos && $next_pos < $l_ref);

        #ignore seq string

        #ignore quality string

        #ignore aux values for now

	$buffer->{POS} = $i;
        return $i; #seems to be a valid bam line
    }

    return;
}

sub __inflate_to_buffer {
    my $self = shift;
    my $buffer = $self->__buffer;

    #remove passed blocks
    while(@{$buffer->{BLOCKS}}){
	my $block = $buffer->{BLOCKS}[0];
	if($buffer->{POS} >= $block->{I_LENGTH} + 10000){ #past current position by 10kb
	    shift @{$buffer->{BLOCKS}}; #remove block
	    $_->{BUFFER_OFFSET} -= $block->{I_LENGTH} foreach(@{$buffer->{BLOCKS}});
	    substr($buffer->{DATA}, 0, $block->{I_LENGTH}, ''); #remove text from buffer
	    $buffer->{POS} -= $block->{I_LENGTH};
	    $buffer->{LENGTH} -= $block->{I_LENGTH};
	}
	else{
	    last;
	}
    }

    #get current block
    my $block = $self->{BLOCK};
    return unless($block->{DATA});

    #inflate block to buffer
    if($block->{I_LENGTH} != 0 || $block->{FILE_OFFSET} == $self->__file->{LENGTH}-28){ #skip empty non-eof blocks
	$block->{BUFFER_OFFSET} = $buffer->{LENGTH};
	$self->__inflate(\ ($block->{DATA}), \ ($buffer->{DATA}));
	$buffer->{LENGTH} += $block->{I_LENGTH};
	push(@{$buffer->{BLOCKS}}, $block); #makes copy of block
    }

    #reset block
    delete($self->{BLOCK});

    return;
}
