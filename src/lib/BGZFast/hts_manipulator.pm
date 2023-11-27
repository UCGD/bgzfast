=head1 NAME

BGZFast::hts_manipulator - Manipulation of sequence realated file structure

=head1 SYNOPSIS

  use BGZFast::hts_manipulator;

  # Jump to start of first feature after given file position
  my $obj = BGZFast::hts_manipulator->new('/path/to/file.gff.gz');
  my ($real, $virtual) = $obj->smart_seek(10_000_000); #pos is real file

  # Get position from index file
  ($real, $virtual) = $obj->query_index($chr, $pos);
  $obj->block_seek($real, $virtual); #move to position

  #read features one after the other to perform analysis
  while(my $f = $obj->next_feature){
    my ($real, $virtual) = $f->offsets();
  }

=head1 DESCRIPTION

BGZFast::hts_manipulator provides a mechanism to easilly navegate and
manipulate high throughput sequence related file formats that have been indexed.
This allows these BGZF compressed files to be shattered and partitioned for
parallel processing and streaming.

=head1 TBI INDEXING

TBI indexing can be performed using the $obj->make_tbi method

=head1 SEE ALSO

L<BGZFast::hts_utility>
L<BGZFast::bgzf_manipulator>

=head1 AUTHOR

Carson Holt E<lt>carson.holt@genetics.utah.eduE<gt>.

Copyright (c) 2017 University of Utah.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head1

For simple access, the following methods are provided:

=cut

#------------------------------------------------------------------------
#----                 BGZFast::hts_manipulator                  ----
#------------------------------------------------------------------------
package BGZFast::hts_manipulator;

use strict;
use warnings;
use Class::Lstruct;
use Fcntl qw(:DEFAULT :seek); #SEEK_SET, SEEK_CUR, SEEK_END
use BGZFast::hts_utility;

use parent qw(BGZFast::bgzf_manipulator);

struct ('BGZFast::header' => [
    DATA              => '&$',
    TAIL              => '&$',
    TAIL_FILE_OFFSET  => '&$',
    TAIL_BLOCK_OFFSET => '&$',
    D_LENGTH          => '&$',
    TEXT              => '&$',
    N_REF             => '&$',
    NAME2SEQID        => '&%',
    NAME2LENGTH       => '&%',
    REF               => '&@',
]);

struct ('BGZFast::feature' => [
    DATA         => '&$',
    LENGTH       => '&$',
    FILE_OFFSET  => '&$',
    BLOCK_OFFSET => '&$',
    REFID        => '&$',
    POS          => '&$',
    END          => '&$',
    NEXT_REFID   => '&$',
    NEXT_POS     => '&$',
]);

struct ('BGZFast::tbi' => [
    PATH       => '&$',
    ABS_PATH   => '&$',
    TYPE       => '&$',
    COL        => '&@',
    META       => '&$',
    SKIP       => '&$',
    LENGTH     => '&$',
    INDICES    => '&@',
    N_NO_COOR  => '&$',
    NAMES      => '&@',
    NAME2TBID  => '&%',
    SEQID2TBID => '&@',
    TBID2SEQID => '&@',
    REF        => '&@'
]);

our $VERSION = '0.01';

#------------------------------------------------------------------------
#--------------       OBJECT INITIALIZATION METHODS        --------------
#------------------------------------------------------------------------
=head2 new

 Title   : new
 Usage   : my $obj = BGZFast::hts_manipulator->new($path, %options);
 Function: Initialize a new file manipulator object.
 Returns : A new BGZFast::hts_manipulator object.
 Args    : A single BGZF file (or regular file if setting -is_bgzf => 0)
           Optional arguments:

 Option        Description                                         Default
 -----------   -----------                                         -------
 -is_bgzf      Indicates whether a file is in BGZF format          1
 -cpus         How many CPUs to use on paralelizable operations    1
 -idx_type     Predefined index type (SAM, VCF, GENERIC, ZERO)     GENERIC
 -idx_col      Seq, start, and end columns to use for indexing     [1,4,5]
 -idx_meta     Meta data character to use when parsing file        #
 -idx_skip     Count of lines to always skip for header            0

=cut
sub new {
    return shift->SUPER::new(@_);
}
#------------------------------------------------------------------------
sub __init {
    my ($self, @args) = @_;

    #initialize base class parameters (returns remaining params)
    my $param = $self->SUPER::__init(@args);

    #class specific parameters
    if(defined(my $arg = delete($param->{IDX_TYPE}))){
        if(uc($arg) eq 'GENERIC' || $arg eq 0){ #generic
	    $self->{IDX_TYPE} = 0;
            $self->{IDX_COL} = [1,4,5];
            $self->{IDX_META} = '#';
            $self->{IDX_SKIP} = 0;
	    $self->{IDX_COLMAX} = (sort {$b <=> $a} @{$self->{IDX_COLS}})[0];
        }
        elsif(uc($arg) eq 'SAM' || $arg eq 1){ #sam
	    $self->{IDX_TYPE} = 1;
            $self->{IDX_COL} = [3,4,0];
            $self->{IDX_META} = '@';
            $self->{IDX_SKIP} = 0;
	    $self->{IDX_COLMAX} = 6;
        }
        elsif(uc($arg) eq 'VCF' || $arg eq 2){ #vcf
	    $self->{IDX_TYPE} = 2;
            $self->{IDX_COL} = [1,2,0];
            $self->{IDX_META} = '#';
            $self->{IDX_SKIP} = 0;
	    $self->{IDX_COLMAX} = 8;
        }
        elsif(uc($arg) eq 'ZERO' || $arg eq 65536){ #generic-zero
	    $self->{IDX_TYPE} = 65536;
            $self->{IDX_COL} = [1,2,3];
            $self->{IDX_META} = '#';
            $self->{IDX_SKIP} = 0;
	    $self->{IDX_COLMAX} = (sort {$b <=> $a} @{$self->{IDX_COLS}})[0];
        }
	else{
	    die "ERROR: Invalid value for IDX_TYPE\n";
	}

	#get index parameters for generic types
	if($self->{IDX_TYPE} == 0 || $self->{IDX_TYPE} == 65536){
	    if(defined(my $arg = delete($param->{IDX_COL}))){
		die "ERROR: Invalid value for IDX_COL\n" if(ref($arg) ne 'ARRAY' || @$arg != 3);
		$self->{IDX_COL} = $arg;
		$self->{IDX_COLMAX} = (sort {$b <=> $a} @{$self->{IDX_COLS}})[0];
	    }
	    if(defined(my $arg = delete($param->{IDX_META}))){
		die "ERROR: Invalid value for IDX_META\n" if(length($arg) != 1);
		$self->{IDX_META} = $arg;
	    }
	    if(defined(my $arg = delete($param->{IDX_SKIP}))){
		die "ERROR: Invalid value for IDX_SKIP\n" if($arg < 0);
		$self->{IDX_SKIP} = $arg;
	    }
	}
	else{ #ignore values on predefined types
	    delete($param->{IDX_COL});
	    delete($param->{IDX_META});
	    delete($param->{IDX_SKIP});
	}
    }

    return $param;
}
#------------------------------------------------------------------------
#used by copy constructor
sub __args {
    my $self = $_[0];

    my @args = $self->SUPER::__args();
    push(@args, IDX_TYPE => $self->{IDX_TYPE});
    push(@args, IDX_COL  => [@{$self->{IDX_COL}}]);
    push(@args, IDX_META => $self->{IDX_META});
    push(@args, IDX_SKIP => $self->{IDX_SKIP});

    return @args;
}
#------------------------------------------------------------------------
#-------------- PLACEHOLDER AND OBJECT METHODS TO OVERRIDE --------------
#------------------------------------------------------------------------
#add feature as a position sensitive value
sub __reset_seek_sensitive {
    my $self = $_[0];
    $self->SUPER::__reset_seek_sensitive();
    $self->{FEATURE} = undef;
}

#------------------------------------------------------------------------
#--------------           HEADER RELATED METHODS           --------------
#------------------------------------------------------------------------
sub __header_struct {
    return BGZFast::header->new(DATA              => '',
				TAIL              => '',
				TAIL_FILE_OFFSET  => undef,
				TAIL_BLOCK_OFFSET => undef,
				D_LENGTH          => 0,
				TEXT              => '',
				N_REF             => undef,
				NAME2SEQID        => {},
				NAME2LENGTH       => {},
				REF               => []);
}
#------------------------------------------------------------------------
sub __header {
    my $self = $_[0];
    $self->{HEADER} = $self->__parse_header if(!$self->{HEADER});
    return $self->{HEADER};
}
#------------------------------------------------------------------------
=head2 header

 Title   : header
 Usage   : my $header = $obj->header();
 Function: Get header from file.
 Returns : A BGZFast::header object.
 Args    : None

=cut
sub header {
    my $self = $_[0];
    return $self->__header->text;
}
#------------------------------------------------------------------------
#parse header from file
sub __parse_header {
    my $self = $_[0];
    my $header = $self->__header_struct;

    #parse reference info out of header
    $self->__fill_header_data($header); #separates header from features
    $self->__fill_header_ref($header); #parses header to fill ref info
    
    return $header;
}
#------------------------------------------------------------------------
#separates header from feature data
sub __fill_header_data {
    my $self = $_[0];
    my $header = $_[1];
    my $IN = $self->__handle;

    #backup current position then move to start of file
    my $c_real_pos = sysseek($IN, 0, SEEK_CUR);
    sysseek($IN, 0, SEEK_SET);

    #read header
    my @h_blocks;
    my $io_buffer = '';
    my $io_offset = 0;
    my $h_buffer = '';
    my $h_offset = 0;

    #skip lines in header
    my $skip = $self->{IDX_SKIP};
    while($skip){
	#find end of line
	my $pos = index($h_buffer, "\n", $h_offset)+1;

	#read in more data if needed
	if(!$pos){
	    #get block
	    my $block = $self->__block_from_ref(\$io_buffer, $io_offset);
	    if(!$block){
		my $stat = $self->__read($IN, $io_buffer, $self->max_bsize, length($io_buffer));
		next if($stat > 0);
	    }
	    die "ERROR: Malformed file. Missing lines.\n" if(!$block);
	    $io_offset += $block->d_length;
	    push(@h_blocks, $block);
	    
	    #grow buffer
	    $block->buffer_offset = length($h_buffer);
	    if($self->{IS_BGZF}){
		$self->__inflate(\ ($block->data), \$h_buffer);
	    }
	    else{
		$h_buffer .= $block->data;
		$block->data = '';
	    }
	    next;
	}

	#skip line
        $skip--;
        $h_offset = $pos;
    }

    #read comment lines
    my $meta = $self->{IDX_META};
    while(defined($meta)){
	#get leading charcter of next line
	my $lead = substr($h_buffer, $h_offset, 1);
	my $pos = index($h_buffer, "\n", $h_offset)+1;

	#read in more data if needed
	if(!length($lead) || ($lead eq $meta && !$pos)){
	    #get block
	    my $block = $self->__block_from_ref(\$io_buffer, $io_offset);
	    if(!$block){
		my $stat = $self->__read($IN, $io_buffer, $self->max_bsize, length($io_buffer));
		next if($stat > 0);
	    }
	    last if(!$block && !length($lead)); #end of file
	    die "ERROR: Malformed file. Missing endline.\n" if(!$block);
	    $io_offset += $block->d_length;
	    push(@h_blocks, $block);

	    #grow buffer
	    $block->buffer_offset = length($h_buffer);
	    if($self->{IS_BGZF}){
		$self->__inflate(\ ($block->data), \$h_buffer);
	    }
	    else{
		$h_buffer .= $block->data;
		$block->data = '';
	    }
	    next;
	}

	#skip line
	last if($lead ne $meta);
	$h_offset = $pos;
    }

    #separate header from feature tail
    my $h_data = substr($h_buffer, 0, $h_offset);
    $header->text = $h_data;

    #identify real and virtual offset for end of header (start of features)
    if(!$self->{IS_BGZF}){ #uncompressed files
        $header->tail_file_offset = length($h_data); #real offset
        $header->tail_block_offset = 0; #virtual offset
        $header->d_length = length($h_data); #trim
        $h_buffer = ''; #no tail
	$header->data = $h_data;
	$header->tail = $h_buffer;
    }
    else{
	my ($block) = grep {$h_offset < $_->buffer_offset + $_->i_length} @h_blocks;
	if($block){
	    my $tb_off = $h_offset-$block->buffer_offset;
	    $header->tail_file_offset = $block->file_offset; #real offset
	    $header->tail_block_offset = $tb_off; #virtual offset
	    $h_buffer = substr($h_buffer, $h_offset, $block->i_length-$tb_off);
	}
	else{
	    $header->tail_file_offset = $io_offset; #real offset
	    $header->tail_block_offset = 0; #virtual offset
	    $h_buffer = ''; #no tail
	}
	$self->zip($h_data, $header->data);
	$self->zip($h_buffer, $header->tail);
    }
    
    #restore file positions
    sysseek($IN, $c_real_pos, SEEK_SET);
    
    return $header;
}
#------------------------------------------------------------------------
#parses header reference text (should be overridden by child module)
sub __fill_header_ref {
    my $self = $_[0];
    my $header = $_[1];

    my $name2seqid = $header->name2seqid;
    my $name2length = $header->name2length;

    #add unmapped
    $name2seqid->{'*'} = -1; #unmapped chr
    $name2length->{'*'} = 0; #unmapped chr

    #Example
    #my $text = $header->text; #header raw text
    #...
    #$header->n_ref = $n_ref;  #total reference seq IDs
    #my $ref = $header->ref;
    #for(my $i = 0; $i < $n_ref; $i++){
    #	$ref->[$i]{NAME}   = $name; #name of seq
    #	$ref->[$i]{LENGTH} = $l_ref; #length of seq
    #	$name2seqid->{$name}  = $i; #name2seqid translation index
    #	$name2length->{$name} = $l_ref; #name2length translation index
    #}

    return $header;
}

#------------------------------------------------------------------------
#--------------           FEATURE RELATED METHODS          --------------
#------------------------------------------------------------------------
sub __feature_struct {
    return BGZFast::feature->new(DATA         => '',
				 LENGTH       => 0,
				 FILE_OFFSET  => undef,
				 BLOCK_OFFSET => undef,
				 REFID        => undef,
				 POS          => undef,
				 END          => undef,
				 NEXT_REFID   => undef,
				 NEXT_POS     => undef);
}
#------------------------------------------------------------------------
sub __feature {
    my $self = $_[0];
    return $self->{FEATURE};
}
#------------------------------------------------------------------------
=head2 feature_stat

 Title   : feature_stat
 Usage   : ($length,
            $refid,
            $pos,
            $end,
            $next_refid,
            $next_pos) = $obj->feature_stat(\$string);
 Function: Get basic feature info from string
 Returns : Feature length, ref ID, beginning pos, end pos, next ref ID for
           split features, next pos
           (returns undef for partial features)
 Args    : String or string reference
           Optional offset in string

=cut
sub feature_stat {
    my $self = $_[0];
    my $ref  = (ref($_[1])) ? $_[1] : \ ($_[1]);
    my $boff = $_[2] || 0;
    my $type = $self->{IDX_TYPE};   #idx_type
    my $c    = $self->{IDX_COL};    #idx_col
    my $meta = $self->{IDX_META};
    my $max  = $self->{IDX_COLMAX};

    #set info for predefined types
    if($type == 1){ #sam format gets end from cigar column
        $max = 6;
	$meta = '@';
    }
    elsif($type == 2){ #vcf format gets end from ref and info column
        $max = 8;
	$meta = '#';
    }
    return if(substr($$ref, $boff, 1) eq $meta); #comment line

    #find end offset to stop search
    my $eoff = index($$ref, "\n", $boff)+1;
    return undef if(!$eoff); #not a full feature
    my $flen = $eoff-$boff;
    $eoff--; #set end offset before endline

    #parse tabs
    my $i = 1;
    my @t = ($boff); #tab positions (offset, length)
    while(1){
        $t[$i*2] = index($$ref, "\t", $t[$i*2-2])+1; #set offset
        $t[$i*2-1] = $t[$i*2] - $t[$i*2-2]-1; #set prev length
        if(!$t[$i*2] || $t[$i*2] > $eoff){ #fix if past end
            shift(@t); #drop partial
            $t[$i*2-1] = $eoff - $t[$i*2-2];
            last;
        }
        elsif(defined($max) && $i == $max){ #drop partial
            shift(@t);
	    last;
        }
        $i++;
    }

    #get stats for type
    my ($refid, $pos, $end, $next_refid, $next_pos);
    if($type == 1){ #SAM
	$refid = $self->name2seqid(substr($$ref, $t[2*2], $t[2*2+1])); #col 3

	#get position
	if($refid == -1){ #unmapped
	    $pos = -1;
	    $end = -1;
	}
	else{ #mapped
	    $pos = substr($$ref, $t[3*2], $t[3*2-1])-1; #col 4

	    #getting end from cigar string
	    my $cigar = substr($$ref, $t[5*2], $t[5*2+1]); #col 6
	    if($cigar ne '*'){
		my $op;
		my $off = 0;
		$end = $pos;
		for(my $i = 0; $i < length($cigar); $i++){
		    $op = substr($cigar, $i, 1);
		    if(ord($op) > 57){ #not a number
			$op = uc($op); #make upper case
			$end += substr($cigar, $off, $i-$off)
			    if($op eq 'M' || $op eq 'D' || $op eq 'N');
			$off = $i+1; #set start offset of next number
		    }
		}
	    }
	    else{ #assume length of 1
		$end = $pos+1;
	    }
	}

	#get next read in pair
	my $next_ref = substr($$ref, $t[6*2], $t[6*2+1]); #col 7
	$next_refid = ($next_ref eq '=') ? $refid : $self->name2seqid($next_ref);
	$next_pos = ($next_refid == -1) ? -1 : substr($$ref, $t[7*2], $t[7*2-1])-1; #col 8
    }
    elsif($type == 2){ #VCF
	$refid = $self->name2seqid(substr($$ref, $t[0*2], $t[0*2+1])); #col 1
        $pos = substr($$ref, $t[1*2], $t[1*2-1])-1; #col 2 (-1 to make zero based)
	
	#getting end is a little more complicated
	my $info = substr($$ref, $t[7*2], $t[7*2+1]); #col 8
	if(substr($info, 0, 4) eq 'END='){ #set using info column if END= is present
	    my $l = index($info, ';', 4)-4;
	    $l = length($info)-4 if($l < 0);
	    $end = substr($info, 4, $l);
	}
	elsif(my $o = index($info, ';END=')+1){
	    $o += 4; #set offset just after ';END='
	    my $l = index($info, ';', $o)-$o;
	    $l = length($info)-$o if($l < 0);
	    $end = substr($info, $o, $l);
	}
	elsif(substr($$ref, $t[3*2], 1) eq '<'){ #check for angle bracket <ID>
	    die "ERROR: Support for <ID> in VCF REF column not implemented\n";
	}
	else{ #set based on ref column length
	    $end = $pos+$t[3*2+1]; #col4
	}
	
	#add breakpoints for SVs in VCF
	my $alt = substr($$ref, $t[4*2], $t[4*2+1]); #col 5
	if(index($alt, '[', 0)+1 || index($alt, ']', 0)+1){ #SV symbols used
	    $alt =~ /[\[\]]\s*([^\[\]\:\s]+)\s*\:\s*([^\[\]\:\s]+)\s*[\[\]]/; #slow
	    $next_pos = $2-1; #make zero based
	    $next_refid = $self->name2seqid($1);
	    if(index($alt, ',', 0)+1){ #too complex
		warn "WARNING: Feature is too complex:\n".substr($$ref, $boff, $eoff-$boff)."\n";
		$next_pos   = -2;
		$next_refid = -2;
	    }
	}
    }
    else{ #Generic and Zero-based
	$refid = $self->name2seqid(substr($$ref, $t[$c->[0]*2-2], $t[$c->[0]*2-1]));
	if($refid == -1){
	    $pos = -1;
	    $end = -1;
	}
	else{
	    $pos = substr($$ref, $t[$c->[1]*2-2], $t[$c->[1]*2-1]);
	    $pos-- unless($type & 65536); #make zero based
	    $end = ($c->[2] && $c->[2] != $c->[1]) ?
		substr($$ref, $t[$c->[2]*2-2], $t[$c->[2]*2-1])-1 : $pos+1;
	}
	$next_pos   = -1;
	$next_refid = -1;
    }

    return ($flen, $refid, $pos, $end, $next_refid, $next_pos);
}
#------------------------------------------------------------------------
=head2 next_feature

 Title   : next_feature
 Usage   : $f = $obj->next_feature();
 Function: Get next feature from buffer
 Returns : Feature object (undef if end of file)
 Args    : None

=cut
sub next_feature {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;

    $self->__read_feature();

    return $self->__feature;
}

#------------------------------------------------------------------------
#walks across uncompressed block looking for feature start
sub __adjust_to_next_feature {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;
    my $meta = $self->{IDX_META};
    my $b_ref = \ ($buffer->data);

    #already at start of feature
    if($buffer->feature_ok){
	my $pos = $buffer->pos;
	my $size = $buffer->length;
	if(length($meta) && ($pos == $size || $buffer->lengthsubstr($$b_ref, $pos, 1) eq $meta)){
	    $buffer->feature_ok = 0;
	}
	else{
	    return $buffer->pos;
	}
    }

    #move past header if necessary
    my ($poff, $voff) = $self->block_tell();
    my ($hpoff, $hvoff) = ($self->__header->tail_file_offset, $self->__header->tail_block_offset);
    $self->block_seek($hpoff, $hvoff) if($poff < $hpoff || ($poff == $hpoff && $voff < $hvoff));

    #get first encountered feature offset
    while(1){
	my $pos = $buffer->pos;
	my $size = $buffer->length;
	if($pos == $size){ #need to grow buffer
	    last unless($self->grow_vbuffer());
	    next;
	}
	$pos = $self->__find_feature_offset($b_ref, $pos, $buffer->upstream);
	if($pos == -1){ #keep looking from end of buffer
	    $buffer->pos = $size;
	    next;
	}

	#check for and skip past meta lines
	if(length($meta) && substr($$b_ref, $pos, 1) eq $meta){
	    $buffer->pos = $pos + 1;
	    next;
	}

	$buffer->pos = $pos;
	$buffer->feature_ok = 1;
	last;
    }

    return $buffer->pos;
}
#------------------------------------------------------------------------
sub __find_feature_offset {
    my $self = $_[0];
    my $offset = $_[2] || 0;
    my $buffer = $self->__vbuffer;
    my $b_ref = \ ($buffer->data);
    my $size = $buffer->length;
    my $meta = $self->{IDX_META};
    
    $offset = 0 if($offset < 0);

    #find endline in data
    while(1){
	my $pos;
	my $upstream = ($offset == 0) ? $buffer->upstream : substr($$b_ref, $offset-1, 1);
	if($upstream eq "\n"){
	    $pos = $offset;
	}
	else{
	    $pos = index($$b_ref, "\n", $offset)+1; #position after endline
	    return -1 if($pos == 0);
	}

	#skip meta lines
	if(length($meta)){
	    if($pos == $size){
		return -1; #can't be certain without next character
	    }
	    elsif(substr($$b_ref, $pos, 1) eq $meta){ #is meta
		$offset = $pos+1;
		next;
	    }
	}

	return $pos;
    }
}

#------------------------------------------------------------------------
#assumes buffer position is already set to start of feature
sub __read_feature {
    my $self = $_[0];
    my $IN = $self->__handle;

    #make feature
    my $b_ref = $self->vbuffer;
    my $b_pos = $self->__adjust_to_next_feature();
    my $feat = $self->__feature_from_ref($b_ref, $b_pos); #make feature
    while(!$feat){
	last if(!$self->grow_vbuffer());
	$b_pos = $self->__adjust_to_next_feature();
	$feat = $self->__feature_from_ref($b_ref, $b_pos);
    }
    $self->{FEATURE} = $feat;
    return 0 if(!$feat);

    #set values
    ($feat->file_offset, $feat->block_offset) = $self->block_tell();
    $self->vbuffer_offset($feat->length, SEEK_CUR); #move past current feature

    return 1;
}
#------------------------------------------------------------------------
#assumes file position is set to start of block
sub __feature_from_ref {
    my $self = $_[0];
    my $ref = ref($_[1]) ? $_[1] : \$_[1];
    my $offset = $_[2] || 0;

    return undef if(length($$ref) == 0);

    #fill in feature
    my $feat = $self->__feature_struct();
    ($feat->length,
     $feat->refid,
     $feat->pos,
     $feat->end,
     $feat->next_refid,
     $feat->next_pos) = $self->feature_stat($ref, $offset);

    return undef if(!defined($feat->length));
    $feat->data = substr($$ref, $offset, $feat->length);

    return $feat;
}

#------------------------------------------------------------------------
#--------------            INDEX RELATED METHODS           --------------
#------------------------------------------------------------------------
sub idx {return tbi(@_)} #sugar
sub tbi {
    my $self = $_[0];
    
    if(!$self->{TBI}){
	$self->{TBI} = (-f $self->tbi_file) ? $self->parse_tbi() : $self->make_tbi();
    }
    return $self->{TBI};
}
#------------------------------------------------------------------------
sub __tbi {
    my $self = $_[0];
    return $self->{TBI};
}
#------------------------------------------------------------------------
sub __tbi_struct { #overrides parent function
    return BGZFast::tbi->new(PATH       => undef,
			     ABS_PATH   => undef,
			     TYPE       => 0,       #0:generic; 1:SAM; 2:VCF, 65536:zero-based
			     COL        => [1,4,5], #column for seq, beg, end
                             META       => '#',     #comment character
			     SKIP       => 0,       #number of lines to skip at beginning
			     LENGTH     => 0,
			     INDICES    => [],
			     N_NO_COOR  => 0,
			     NAMES      => [],
			     NAME2TBID  => {},
			     SEQID2TBID => [],
			     TBID2SEQID => [],
			     REF        => []);
}
#------------------------------------------------------------------------
sub index_file {return(tbi_file(@_))} #sugar
sub tbi_file { #overrides parent function
    my $self = $_[0];
    my $file = $_[1];

    if($file){
        return $file.".tbi";
    }
    elsif($self->{TBI} && defined($self->{TBI}->path)){
        return $self->{TBI}->path;
    }
    return $self->file.".tbi";
}
#------------------------------------------------------------------------
sub make_index {return make_tbi(@_)}
sub make_tbi {
    my $self = $_[0];

    my $file = $self->file;
    my $tbi_file = $self->tbi_file;
    if(my $exe = File::Which::which('tabix')){
	#build command for index type => 0: generic (1 based); 1: SAM; 2: VCF, 65536: generic (zero based)
	my @cmd = ($exe, '-f');
	if($self->{IDX_TYPE} == 1){
	    push(@cmd, '-p' => 'sam');
	}
	elsif($self->{IDX_TYPE} == 2){
	    push(@cmd, '-p' => 'vcf');
	}
	else{
	    push(@cmd,
		 '-s' => $self->{IDX_COL}[0], #seq col
		 '-b' => $self->{IDX_COL}[1], #begin col
		 '-e' => $self->{IDX_COL}[2], #end col
		 '-c' => $self->{IDX_META}, #comment/meta character
		 '-S' => $self->{IDX_SKIP}); #fixed numer of header lines to skip
	    push(@cmd, '-0') if($self->{IDX_TYPE} == 65536);
	}
	push(@cmd, $file);
	
	#run indexer
	my $stat = system(@cmd);
	unlink($tbi_file) if($stat); #on failure
    }
    die "ERROR: ".__PACKAGE__." could not generate TBI index files\n" if(!-f $tbi_file);

    return $self->parse_tbi($tbi_file);
}
#------------------------------------------------------------------------
sub parse_tbi {
    my $self = $_[0];
    my $tbi_file = $_[1] || $self->tbi_file;
    
    #read in entire file
    my $raw; #raw data
    my $size = (stat($tbi_file))[7];
    sysopen(my $IN, $tbi_file, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $tbi_file: $!";
    $self->__read($IN, $raw, $size);
    close($IN);

    #inflate tbi (it's BGZF compressed)
    my $data;
    while(length($raw)){
        #read BGZF header
        my $b_data = substr($raw, 0, 18);
        die "ERROR: Failure to read BGZF string\n" if(length($b_data) != 18);

        #validate
        my ($id1,$id2,$cm,$flg,$mtime,$xfl,$os,
            $xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $b_data);
        die "ERROR: Does not appear to be a BGZF file\n"
            if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
               $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);

        #read remainder compression block and footer
        my $c_data_size = $bsize-$xlen-19; #compression block size
        $b_data .= substr($raw, 0, $c_data_size + 8, ''); #the +8 is for the footer
        die "ERROR: Could not read compression block\n"
            if(length($b_data) != 18 + $c_data_size + 8);

        #inflate and concatenate to output
        $self->__inflate(\$b_data, \$data);
    }

    my $tbi = $self->{TBI} = $self->__string2tbi($data);
    $tbi->path = $tbi_file;
    $tbi->abs_path = Cwd::abs_path($tbi_file);
    $tbi->length = $size;
    
    return $tbi;
}
#------------------------------------------------------------------------
sub __string2tbi {
    my $self = $_[0];
    
    #check magic string
    my $str_off = 0; #offset in substring
    if(substr($_[1], $str_off, 4) ne "TBI\1"){
        die "ERROR: File does not appear to be a TBI index file\n"
    }
    $str_off += 4;

    #get number of reference sequences
    my $n_ref = unpack('l<', substr($_[1], $str_off, 4));
    $str_off += 4;

    #get formatting info
    my $format = unpack('l<', substr($_[1], $str_off, 4));
    $str_off += 4;
    my ($col_seq, $col_beg, $col_end) = unpack('l<l<l<', substr($_[1], $str_off, 12));
    $str_off += 12;
    my ($meta, $skip) = unpack('l<l<', substr($_[1], $str_off, 8));
    $meta = chr($meta); #convert to character
    $str_off += 8;

    #get ref info
    my $l_nm = unpack('l<', substr($_[1], $str_off, 4));
    $str_off += 4;
    my @names = grep {length($_)} split(/\0/, substr($_[1], $str_off, $l_nm));
    my %name2tbid = map {$names[$_] => $_} (0..$#names);
    $str_off += $l_nm;

    #list of indices (n=n_ref)
    my @tbi; #list representation
    for(my $i = 0; $i < $n_ref; $i++){
        $tbi[$i] = {BINS => [], LINEAR => [], PSEUDO => undef}; #initialize index of each contig

        #get number of distinct bins (for the binning index)
        my $n_bin = unpack('l<', substr($_[1], $str_off, 4));
        $str_off += 4;

	#list of distinct bins (n=n_bin)
	my $tbi_b = $tbi[$i]{BINS};
        for(my $j = 0; $j < $n_bin; $j++){
            #get distinct bin & number of chunks
            my ($bin, $n_chunk) = unpack('L<l<', substr($_[1], $str_off, 8));
            $str_off += 8;

            #fill in pseudo-bin
            if($bin == 37450){
                my $tbi_p = $tbi[$i]{PSEUDO} = []; #initialize chunk list for pseudo bin
                die "ERROR: Corrupt TBI index for pseudo bin\n" if($n_chunk != 2); #always 2

                #get (virtual) file offset of the start and end of placed unmapped reads
                my ($unmapped_beg, $unmapped_end) = unpack('Q<Q<', substr($_[1], $str_off, 16));
                $str_off += 16;

                #offsets are encoded as --> coffset<<16|uoffset
                my $coffset_beg = ($unmapped_beg >> 16);
                my $uoffset_beg = $unmapped_beg & 65535; #mask is (1<<16)-1
                my $coffset_end = ($unmapped_end >> 16);
                my $uoffset_end = $unmapped_end & 65535; #mask is (1<<16)-1

                push(@{$tbi_p}, [$coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end]);

                #get number of mapped/unmapped read-segments for this reference
                my ($n_mapped, $n_unmapped) = unpack('Q<Q<', substr($_[1], $str_off, 16));
                $str_off += 16;

                push(@{$tbi_p}, [$n_mapped, $n_unmapped]);

                next;
            }

            #list of chunks (n=n_chunk)
            $tbi_b->[$bin] = []; #initialize chunk list for each bin
            for(my $k = 0; $k < $n_chunk; $k++){
                #get (virtual) file offset of the start and end of the chunk
                my ($chunk_beg, $chunk_end) = unpack('Q<Q<', substr($_[1], $str_off, 16));
                $str_off += 16;

                #chunk offsets are encoded as --> coffset<<16|uoffset
                my $coffset_beg = ($chunk_beg >> 16);
                my $uoffset_beg = $chunk_beg & 65535; #mask is (1<<16)-1
                my $coffset_end = ($chunk_end >> 16);
                my $uoffset_end = $chunk_end & 65535; #mask is (1<<16)-1

                push(@{$tbi_b->[$bin]}, [$coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end]);
            }
        }

        #get number of 16kbp intervals (for the linear index)
        my $n_intv = unpack('l<', substr($_[1], $str_off, 4));
        $str_off += 4;

        #list of intervals (n=n_intv)
        my $tbi_l = $tbi[$i]{LINEAR};
        for(my $j = 0; $j < $n_intv; $j++){
            #get (virtual) file offset of the first feature in the interval
            my $ioffset = unpack('Q<', substr($_[1], $str_off, 8));
            $str_off += 8;

            #interval offsets are encoded as --> coffset<<16|uoffset
            my $i_coffset = ($ioffset >> 16);
            my $i_uoffset = $ioffset & 65535; #mask is (1<<16)-1

            $tbi_l->[$j] = [$i_coffset, $i_uoffset];
        }
    }

    #get number of unplaced unmapped reads (RNAME *)
    my $n_no_coor = unpack('Q<', substr($_[1], $str_off, 8));
    $str_off += 8;

    #hash representation
    my $tbi = $self->__tbi_struct();
    $tbi->type = $format;
    $tbi->col = [$col_seq, $col_beg, $col_end];
    $tbi->meta = $meta;
    $tbi->skip = $skip;
    $tbi->indices = \@tbi;
    $tbi->n_no_coor = $n_no_coor;
    $tbi->names = \@names;
    $tbi->name2tbid = \%name2tbid;

    #add usful relation info from file header
    my $ref = $self->__header->ref;
    $tbi->ref = $ref;
    $tbi->seqid2tbid = [map {$name2tbid{$self->seqid2name($_)}} (0..$#$ref)];
    $tbi->tbid2seqid = [map {$self->name2seqid($names[$_])} (0..$#names)];
    
    return $tbi;
}
#------------------------------------------------------------------------
sub write_tbi {
    my $self = $_[0];
    my $tbi_file = $_[1] || $self->tbi_file();

    my $tbi_str = $self->__tbi2string;    

    #compress index
    my $data = '';
    $self->__deflate(\ (substr($tbi_str, 0, 65280, '')), \$data, 6) while(length($tbi_str));
    $data .= $self->eof_block;

    sysopen(my $OUT, $tbi_file, O_RDWR|O_CREAT|O_TRUNC)
	or die "ERROR: Could not open $tbi_file for writing: $!";
    $self->__write($OUT, $data) or die "ERROR: Failed writing to $tbi_file: $!";
    close($OUT);
    
    #update index file info
    my $tbi = $self->__tbi();
    $tbi->path = $tbi_file;
    $tbi->abs_path = Cwd::abs_path($tbi_file);
    $tbi->length = (stat($tbi_file))[7];

    return 1;
}
#------------------------------------------------------------------------
sub __tbi2string {
    my $self = $_[0];
    my $tbi = $_[1] || $self->tbi;

    #add magic
    my $data = "TBI\1";

    #number of reference sequences
    my $n_ref = scalar(@{$tbi->names});
    $data .= pack('l<', $n_ref);

    #formatting info
    $data .= pack('l<', $tbi->type);
    $data .= pack('l<l<l<', @{$tbi->col});
    $data .= pack('l<', ord($tbi->meta));
    $data .= pack('l<', $tbi->skip);

    #ref info
    my $names = join("\0", @{$tbi->names})."\0";
    my $l_nm = length($names);
    $data .= pack('l<', $l_nm);
    $data .= $names;

    #list of indices (n=n_ref)
    for(my $i = 0; $i < $n_ref; $i++){
        #number of distinct bins (for the binning index)
        my $tbi_b = $tbi->indices($i)->{BINS}; #real bins
        my $tbi_p = $tbi->indices($i)->{PSEUDO}; #pseudo bin
        my @bins = grep {$tbi_b->[$_]} (0..$#$tbi_b);
        push(@bins, 37450) if($tbi_p); #optional pseudo bin
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
                my ($coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end) = @{$tbi_p->[0]};
                my $unmapped_beg = ($coffset_beg << 16)|$uoffset_beg;
                my $unmapped_end = ($coffset_end << 16)|$uoffset_end;

                #(virtual) file offset of the start and end of placed unmapped reads
                $data .= pack('Q<Q<', $unmapped_beg, $unmapped_end);

                #number of mapped/unmapped read-segments for this reference
                my ($n_mapped,$n_unmapped) = @{$tbi_p->[1]};
                $data .= pack('Q<Q<', $n_mapped, $n_unmapped);

                next;
            }

            #list of chunks (n=n_chunk)
            my $n_chunk = scalar(@{$tbi_b->[$bin]});
            $data .= pack('L<l<', $bin, $n_chunk);
            for(my $k = 0; $k < $n_chunk; $k++){
                #chunk offsets are encoded as --> coffset<<16|uoffset
                my ($coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end) = @{$tbi_b->[$bin][$k]};
                my $chunk_beg = ($coffset_beg << 16)|$uoffset_beg;
                my $chunk_end = ($coffset_end << 16)|$uoffset_end;

                #(virtual) file offset of the start and end of the chunk
                $data .= pack('Q<Q<', $chunk_beg, $chunk_end);
            }
        }

        #number of 16kbp intervals (for the linear index)
        my $tbi_l = $tbi->indices($i)->{LINEAR};
        my $n_intv = scalar(@$tbi_l);
        $data .= pack('l<', $n_intv);

        #list of intervals (n=n_intv)
        for(my $j = 0; $j < $n_intv; $j++){
            #interval offsets are encoded as --> coffset<<16|uoffset
            my ($i_coffset, $i_uoffset) = @{$tbi_l->[$j]};
            my $ioffset = ($i_coffset << 16)|$i_uoffset;

            #(virtual) file offset of the first feature in the interval
            $data .= pack('Q<', $ioffset);
        }
    }

    #number of unplaced unmapped reads (RNAME *)
    my $n_no_coor = $tbi->n_no_coor;
    $data .= pack('Q<', $n_no_coor) if(defined($n_no_coor));

    return $data;
}
#------------------------------------------------------------------------
sub seqid2name {
    my $self = $_[0];
    my $id = $_[1];
    my $header = $self->__header;

    if(! $header->ref($id)){
	if($self->{TBI}){
	    my $tbi = $self->tbi;
	    if(!defined($tbi->names($id))){
		warn "WARNING: No contig exists with seqID $id\n";
		return undef;
	    }
	    
	    return $tbi->names($id);
	}
	else{
	    warn "WARNING: No contig exists with seqID $id\n";
	    return undef;
	}
    }

    return $header->ref($id)->{NAME};
}
#------------------------------------------------------------------------
sub seqid2seqlength {
    my $self = $_[0];
    my $id = $_[1];
    my $header = $self->__header;

    if(! $header->ref($id)){
	warn "WARNING: No contig exists with seqID $id\n";
	return undef;
    }

    return $header->ref($id)->{LENGTH};
}
#------------------------------------------------------------------------
sub name2seqid {
    my $self = $_[0];
    my $name = $_[1];
    my $header = $self->__header;

    if(!defined($header->name2seqid($name))){
	if($self->{TBI} && !$self->{TBI}{__HEADERFILL}){ #fill from index rather than header
	    my $tbi = $self->tbi;
	    my $names = $tbi->names;
	    for(my $i = 0; $i < @$names; $i++){
		$header->name2seqid($names->[$i]) ||= $i;
	    }
	    $self->{TBI}{__HEADERFILL} = 1; #flag to only use TBI index to fill header once

	    if(!defined($header->name2seqid($name))){
		warn "WARNING: No contig exists with name $name\n";
		return undef;
	    }
	    
	    return $header->name2seqid($name);
	}
	else{
	    warn "WARNING: No contig exists with name $name\n";
	    return undef;
	}
    }

    return $header->name2seqid($name);
}
#------------------------------------------------------------------------
sub name2seqlength {
    my $self = $_[0];
    my $name = $_[1];
    my $header = $self->__header;

    if(! defined($header->name2length($name))){
	warn "WARNING: No contig exists with name $name\n";
	return undef;
    }

    return $header->name2length($name);
}
#------------------------------------------------------------------------
#returns position immediatly prior to unmapped data
sub unmapped_offsets {
    my $self = $_[0];
    my $tbi    = $self->tbi;

    #scan for last mapped bin in index
    my $offset;
    my $v_offset;
    for(my $i = $#{$tbi->indices}; $i >= 0; $i--){
	#use pseudo index first
        if($tbi->indices($i)->{PSEUDO}){
            my $pseudo = $tbi->indices($i)->{PSEUDO}[0];
            ($offset, $v_offset) = ($pseudo->[2], $pseudo->[3]);
        }
        last if(defined($offset));
	
        #use bin index next
        my $bindex = $tbi->indices($i)->{BINS};
        for(my $j = $#$bindex; $j >= 0; $j--){
            my $set = $bindex->[$j];
            next if(!$set || !@$set); #empty bin
	    
            my $l = $set->[-1];
            if(!defined($offset) || $offset < $l->[2] || ($offset == $l->[2] && $v_offset < $l->[3])){
                ($offset, $v_offset) = ($l->[2], $l->[3]);
            }
        }
        last if(defined($offset));
	
        #use linear index next
        my $lindex = $tbi->indices($i)->{LINEAR};
        for(my $j = $#$lindex; $j >= 0; $j--){
            my $set = $lindex->[$j];
            next if($set->[0] == 0 && $set->[1] == 0); #empty bin offset
            ($offset, $v_offset) = ($set->[0], $set->[1]);
            last;
        }
        last if(defined($offset));
    }

    #still nothing? return position immediately after header
    if(!defined($offset)){
        my $h = $self->__header;
        ($offset, $v_offset) = ($h->tail_file_offset, $h->tail_block_offset);
    }

    return ($offset, $v_offset);
}
#------------------------------------------------------------------------
#fast with sacrificed accuracy (assumes feature can only overlap 2 intervals)
sub coverage_estimate {
    my $self = $_[0];
    my $coverage = $_[1] || []; #coverage can be addative
    my $tbi  = $self->tbi;

    #initialize coverage
    if(!@$coverage){
        my $n_ref = $self->__header->n_ref;
        for(my $i = 0; $i < $n_ref; $i++){
            my $len = $self->seqid2seqlength($i);
            my $end_bin = pos2bin($len-1);
            $coverage->[$i] = [map {0} (4681..$end_bin)];
        }
    }

    #use linear index (faster but less accurate than bin index)
    my ($offset2, $v_offset2) = $self->unmapped_offsets();
    for(my $i = $#{$tbi->indices}; $i >= 0; $i--){
        my $j = $tbi->tbid2seqid($i);
        my $lindex = $tbi->indices($i)->{LINEAR};
        for(my $k = $#$lindex; $k >= 0; $k--){
            #Hack for inconsistency in tabix, so terminal positions can't be 0
            if($k < 2 || $#$lindex-2 < $k){
                $coverage->[$j][$k] += 1;
            }

            my $set = $lindex->[$k];
            #next if($set->[0] == 0 && $set->[1] == 0); #empty bin (already set to 0)

            #add coverage
	    my $cov;
            my ($offset1, $v_offset1) = ($set->[0], $set->[1]);
            if($offset1 == $offset2){
		next if($v_offset2 == $v_offset1); #empty bin (already set to 0)
                $cov = (($v_offset2-$v_offset1)>>3)+1; #1/8 (VCFs compress well)
            }
            else{
                $cov = $offset2-$offset1;
            }
            $coverage->[$j][$k] += $cov;

            #assume some coverage leaks into next bin (this hack is why BIN index not needed)
            if($k != $#$lindex){ #leak forward
                $coverage->[$j][$k+1] += ($cov>>8) + 1; #~64 bytes per 16kb (130bp)
            }
            if($k != 0){ #leak backward
                $coverage->[$j][$k-1] += ($cov>>8) + 1; #~64 bytes per 16kb (130bp)
	    }

            ($offset2, $v_offset2) = ($offset1, $v_offset1); #set new boundary
	}
    }

    return $coverage;
}
#------------------------------------------------------------------------
#--------------            DATA SEEK/TELL METHODS          --------------
#------------------------------------------------------------------------

#jump to a given byte position in the file (will auto-adjust to start of next read in block)
sub smart_seek {
    my $self   = $_[0];
    my $pos    = $_[1];
    my $wence  = $_[2];
    my $file = $self->__file;
    my $header = $self->__header;

    $pos = $self->__file_seek($pos, $wence);
    $pos = $self->__file_seek($file->length-28) if($self->{IS_BGZF} && $pos >= $file->length-28); #jump back from eof block
    $pos = $self->__file_seek($header->tail_file_offset) if($pos < $header->tail_file_offset);
    $self->__adjust_to_next_block();

    #keep reading blocks until feature is found
    while($self->__read_block()){
	$self->__inflate_to_vbuffer();
	last if(defined($self->__adjust_to_next_feature()));
    }

    return $self->block_tell;
}

#------------------------------------------------------------------------
#seek to an exact position
sub pos_seek {
    my $self = $_[0];
    my $chr = $_[1];
    my $pos = $_[2];
    my $header = $self->__header;
    my $tbi    = $self->tbi;

    my $offset;
    my $v_offset;
    my $seqid = $self->name2seqid($chr);
    if($seqid == -1){
	#return offset of eof_block if seeking end of unmapped bin
	if($pos >= 0){
	    return wantarray ? ($self->__file('LENGTH')-28, 0) : $self->__file('LENGTH')-28;
	}

	#set start to be end of last mapped bin
	for(my $i = $#{$tbi->indices}; $i >= 0; $i--){
	    my $index = $tbi->indices($i)->{LINEAR};
	    for(my $j = $#$index; $j >= 0; $j--){
		my $set = $index->[$j];
		next if($set->[0] == 0 && $set->[1] == 0); #empty bin offset
		$offset = $set->[0];
		$v_offset = $set->[1];
		last;
	    }
	    last if(defined($offset));
        }
    }
    else{
	#get linear bin
	my $lin = pos2bin($pos) - 4681;
	die "ERROR: Can't seek to unmapped positions\n" if($lin == -1); #investigate unmapped bin

	#scan linear index to get offset
	for(my $i = $seqid; $i < $header->n_ref; $i++){
	    my $tbid = $tbi->seqid2tbid($i); #convert to tabix ID
	    next if(!defined($tbid));

	    my $index = $tbi->indices($tbid)->{LINEAR};
	    for(my $j = ($i == $seqid) ? $lin : 0; $j < $#$index; $j++){
		my $set = $index->[$j];
		next if($set->[0] == 0 && $set->[1] == 0); #empty bin offset
                $offset = $set->[0];
                $v_offset = $set->[1];
                last;
	    }
	    last if(defined($offset));
	}
	    
	#still nothing? scan for last mapped bin
	if(!defined($offset)){
	    for(my $i = $#{$tbi->indices}; $i >= 0; $i--){
		my $index = $tbi->indices($i)->{LINEAR};
		for(my $j = $#$index; $j >= 0; $j--){
		    my $set = $index->[$j];
		    next if($set->[0] == 0 && $set->[1] == 0); #empty bin offset
		    $offset = $set->[0];
		    $v_offset = $set->[1];
		    last;
		}
		last if(defined($offset));
	    }
		
	    #return offset before eof_block if there is nothing else
	    if(!defined($offset)){
		return wantarray ? ($self->__file('LENGTH')-28, 0) : $self->__file('LENGTH')-28;
	    }
	}
    }

    #when the index is empty (no mapped reads) use the end of header position
    if(!defined($offset) && !defined($v_offset)){
	($offset, $v_offset) = ($header->tail_file_offset, $header->tail_block_offset);
    }

    #go to desired position and then adjust for feature match
    $self->block_seek($offset, $v_offset);
    while(my $feat = $self->next_feature){
	my $seqid2 = $feat->refid;
	my $pos2 = $feat->pos; #zero based
	next if($seqid2 < $seqid && $seqid2 != -1); #too far in front
	next if($seqid2 > $seqid && $seqid == -1); #too far in front
	next if($seqid2 == $seqid && $pos2 < $pos); #too far in front
	return $self->block_seek($feat->file_offset, $feat->block_offset);
    }

    #return offset of eof_block
    return wantarray ? ($self->__file('LENGTH')-28, 0) : $self->__file('LENGTH')-28;
}

#------------------------------------------------------------------------
#--------------            HTS_UTILITY ALIASES             --------------
#------------------------------------------------------------------------
#general methods
sub pseudo_bin { shift @_; return BGZFast::hts_utility::pseudo_bin(@_); }
sub min_binsize { shift @_; return BGZFast::hts_utility::min_binsize(@_); }
sub max_binsize { shift @_; return BGZFast::hts_utility::max_binsize(@_); }
sub bin2bit { shift @_; return BGZFast::hts_utility::bin2bit(@_); }
sub bin2size { shift @_; return BGZFast::hts_utility::bin2size(@_); }
sub bin2shift { shift @_; return BGZFast::hts_utility::bin2shift(@_); }

#linear index methods
sub pos2lin { shift @_; return BGZFast::hts_utility::pos2lin(@_); }
sub reg2lins { shift @_; return BGZFast::hts_utility::reg2lins(@_); }

#bin index methods
sub pos2bin { shift @_; return BGZFast::hts_utility::pos2bin(@_);}
sub pos2bins { shift @_; return BGZFast::hts_utility::pos2bins(@_);}
sub reg2bin { shift @_; return BGZFast::hts_utility::reg2bin(@_);}
sub reg2bins { shift @_; return BGZFast::hts_utility::reg2bins(@_);}

#useful conversions
sub lin2bin { shift @_; return BGZFast::hts_utility::lin2bin(@_);}
sub bin2lin { shift @_; return BGZFast::hts_utility::bin2lin(@_);}
sub bin2lins { shift @_; return BGZFast::hts_utility::bin2lins(@_);}
sub bin2beg { shift @_; return BGZFast::hts_utility::bin2beg(@_);}
sub bin2end { shift @_; return BGZFast::hts_utility::bin2end(@_);}
sub bin2reg { shift @_; return BGZFast::hts_utility::bin2reg(@_);}
sub bin2parent { shift @_; return BGZFast::hts_utility::bin2parent(@_);}
sub bin2parents { shift @_; return BGZFast::hts_utility::bin2parents(@_);}
sub bin2children { shift @_; return BGZFast::hts_utility::bin2children(@_);}

#------------------------------------------------------------------------
#--------------      HOUSEKEEPING/CONVENIENCE METHODS      --------------
#------------------------------------------------------------------------
sub __get_tabs {
    my $self = $_[0];
    my $ref = (ref($_[1])) ? $_[1] : \ ($_[1]);
    my $boff = $_[2] || 0;
    my $max = $_[3];

    my $eoff = index($$ref, "\n", $boff)+1 || length($$ref)+1;
    $eoff--;

    my $i = 1;
    my @t = ($boff); #tab positions (offset, length)
    while(1){
        $t[$i*2] = index($$ref, "\t", $t[$i*2-2])+1; #set offset
        $t[$i*2-1] = $t[$i*2] - $t[$i*2-2]-1; #set prev length
        if(!$t[$i*2] || $t[$i*2] > $eoff){ #fix if past end
            shift(@t); #drop partial
            $t[$i*2-1] = $eoff - $t[$i*2-2];
            last;
        }
        elsif(defined($max) && $i == $max){ #drop partial
            shift(@t);
	    last;
	}
	$i++;
    }

    return \@t;
}
#------------------------------------------------------------------------
1;
