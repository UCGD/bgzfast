=head1 NAME

BGZFast::bgzf_manipulator - Manipulation of BGZF file structure

=head1 SYNOPSIS

  use BGZFast::bgzf_manipulator;

  # Jump to start of a BGZF block using real file
  my $obj = BGZFast::bgzf_manipulator->new('/path/to/file.gz');
  my ($real, $virtual) = $obj->smart_seek(10_000_000); #pos is real file
  my $block = $obj->current_block();

  # Move to position using virtual file
  ($real, $virtual) = $obj->pos_seek(10_000_000); #pos is in virtual file
  $block = $obj->current_block();

  # Get position from GZI index file
  ($real, $virtual) = $obj->query_gzi(10_000_000);
  $obj->block_seek($real, $virtual); #move to position

  # Scan for desired contents and get offsets overlapping content
  my $buffer = $obj->vbuffer();
  $obj->grow_vbuffer();
  my $off = index($$buffer, "\n", $boff);
  $obj->vbuffer_offset($off);
  ($real, $virtual) = $obj->block_tell();

  #read blocks one after the other and write to file handle
  while(my $block = $obj->next_block){
     print $obj->write_block($OUT, $block);
  }

=head1 DESCRIPTION

BGZFast::bgzf_manipulator provides a mechanism to easilly navegate and
manipulate BGZF file structure as opposed to manipulation of content. This allows
BGZF files to be shattered and partitioned for parallel processing and streaming
of blocks.

=head1 GZI INDEXING

GZI indexing can be performed using the $obj->make_gzi method

=head1 SEE ALSO

L<BGZFast::bgzf_utility>

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
#----                  BGZFast::bgzf_manipulator                     ----
#------------------------------------------------------------------------

package BGZFast::bgzf_manipulator;

use strict;
use warnings;
use Class::Lstruct;
use Fcntl qw(:DEFAULT :seek); #SEEK_SET, SEEK_CUR, SEEK_END
use BGZFast::bgzf_utility;

use parent qw();

use constant MAX_BSIZE => BGZFast::bgzf_utility::MAX_BSIZE;
use constant MAX_ISIZE => BGZFast::bgzf_utility::MAX_ISIZE;
use constant EOF_BLOCK => BGZFast::bgzf_utility::EOF_BLOCK;

struct ('BGZFast::file' => [
    NAME      => '&$',
    PATH      => '&$',
    ABS_PATH  => '&$',
    HANDLE    => '&$',
    LENGTH    => '&$',
    POS       => '&$',
    BLOCK_OK  => '&$',
    IO_BUFFER => '&$',
]);

struct ('BGZFast::block' => [
    DATA          => '&$',
    D_LENGTH      => '&$',
    I_LENGTH      => '&$',
    FILE_OFFSET   => '&$',
    BUFFER_OFFSET => '&$',
]);

struct ('BGZFast::vbuffer' => [
    DATA       => '&$',
    UPSTREAM   => '&$',
    LENGTH     => '&$',
    POS        => '&$',
    BLOCKS     => '&@',
    FEATURE_OK => '&$',
]);

struct ('BGZFast::gzi' => [
    PATH     => '&$',
    ABS_PATH => '&$',
    LENGTH   => '&$',
    INDICES  => '&%',
]);

our $VERSION = '0.01';

#------------------------------------------------------------------------
#--------------       OBJECT INITIALIZATION METHODS        --------------
#------------------------------------------------------------------------
=head2 new

 Title   : new
 Usage   : my $obj = BGZFast::bgzf_manipulator->new($path, %options);
 Function: Initialize a new file manipulator object.
 Returns : A new BGZFast::bgzf_manipulator object.
 Args    : A single BGZF file (or regular file if setting -is_bgzf => 0)
           Optional arguments:

 Option        Description                                         Default
 -----------   -----------                                         -------
 -is_bgzf      Indicates whether a file is in BGZF format          1
 -cpus         How many CPUs to use on paralelizable operations    1

=cut
sub new {
    my ($class, @args) = @_;
    return $class->clone() if(ref($class) && !@args); #copy constructor

    my $self = {};
    bless($self, ref($class) || $class);

    #initialize
    my $param = $self->__init(@args); #returns unrecognized parameters

    #fail on unrecognized parameters
    die "ERROR: Invalid option(s) (".join(', ', keys %$param).") passed to ".__PACKAGE__."::new"
        if(keys %$param);

    return $self;
}
#------------------------------------------------------------------------
sub clone {
    my $self = $_[0];
    my @args = $self->__args();
    my $clone = "$self"->new(@args);
    $clone->{GZI} = $self->{GZI}; #share indices between copies
    return $clone;
}
#------------------------------------------------------------------------
sub __init {
    my ($self, @args) = @_;
    
    #get user supplied parameters
    my $param;
    if (ref($args[0]) eq 'HASH') {
	$param = _cap_hash($args[0]);
    }
    elsif(@args == 2 && ref($args[1]) eq 'HASH'){
	$param = _cap_hash({%{$args[1]}, FILE => $args[0]});
    }
    elsif(@args % 2 == 1){
	push(@args, FILE => shift @args);
	$param = _cap_hash({@args});
    }
    else {
        $param = _cap_hash({@args});
    }
    
    #empty some values just in case called on already initialized object
    if(defined($param->{FILE}) && defined($self->{FILE}) &&
       Cwd::abs_path($param->{FILE}) ne Cwd::abs_path($self->file)
	){
	$self->__reset_vbuffer;
	$self->__reset_seek_sensitive;
	delete($self->{BLOCK});
	delete($self->{GZI});
	delete($self->{FILE});
	delete($self->{IS_BGZF});
    }
    
    #compression setting
    if(defined($param->{IS_BGZF})){
	$self->{IS_BGZF} = delete($param->{IS_BGZF});
    }
    
    #get file parameter
    if(defined(my $arg = delete($param->{FILE}))){
	die "ERROR: File $arg does not exist" if(! -f $arg);
	
        #set file info
        my $file = $self->__file;
        $file->path     = $arg;
        ($file->name)   = $arg =~ /([^\/]+)$/;
        $file->abs_path = Cwd::abs_path($arg);
        $file->pos      = 0;
        $file->length   = (stat($arg))[7];
	$file->io_buffer = '';
	
        #open file
        my $IN = $self->__handle;
	
        #detect BGZF if necessary
        $self->{IS_BGZF} = 1 if($arg =~ /\.gz$/ && !defined($self->{IS_BGZF})); #BGZF extension
        if(!defined($self->{IS_BGZF})){
            my $data = '';
            sysseek($IN, 0, SEEK_SET);
            $self->__read($IN, $data, 16, 0);
	    sysseek($IN, 0, SEEK_SET); #restore
	    $self->{IS_BGZF} = 1 if($self->validate_block_header($data));
        }
	
        #validate EOF
        if($self->{IS_BGZF}){
            my $last = '';
            sysseek($IN, -28, SEEK_END); #last 28 bytes
            $self->__read($IN, $last, 28, 0);
            sysseek($IN, 0, SEEK_SET); #restore
            die "ERROR: The file $arg is missing the EOF marker"
                if($last ne $self->eof_block());
        }
    }
    else{
	die "ERROR: You must provide a file to initialize ".__PACKAGE__;
    }

    #get gzi index file parameter
    if(defined(my $arg = delete($param->{GZI_FILE}))){
	my $gzi = $self->__gzi;
	$gzi->path = $arg;
	$gzi->abs_path = Cwd::abs_path($arg);
    }
    
    return $param;
}
#------------------------------------------------------------------------
sub __args {
    my $self = $_[0];

    my @args;
    push(@args, FILE => $self->file);
    push(@args, CPUS => $self->{CPUS});
    push(@args, IS_BGZF => $self->{IS_BGZF});

    return @args;
}
#------------------------------------------------------------------------
#-------------- PLACEHOLDER AND OBJECT METHODS TO OVERRIDE --------------
#------------------------------------------------------------------------

#placeholder for values that need to be cleared on seek operations
sub __reset_seek_sensitive {
    return;
}

#------------------------------------------------------------------------
#--------------            FILE RELATED METHODS            --------------
#------------------------------------------------------------------------
sub __file_struct {
    return BGZFast::file->new(NAME     => undef,
			      PATH     => undef,
			      ABS_PATH => undef,
			      HANDLE   => undef,
			      LENGTH   => 0,
			      POS      => 0,
			      BLOCK_OK => undef,
			      IO_BUFFER  => '');
}
#------------------------------------------------------------------------
sub __file {
    my $self = $_[0];
    $self->{FILE} ||= $self->__file_struct();
    return $self->{FILE};
}
#------------------------------------------------------------------------
=head2 file

 Title   : file
 Usage   : my $path = $obj->file;
 Function: Get the path to the BGZF file encapsulated by this object.
 Returns : String
 Args    : None

=cut

sub file {
    return $_[0]->__file->path;
}
#------------------------------------------------------------------------
sub __handle {
    my $self = $_[0];
    my $file = $self->__file;
    $self->__open_handle if(!$file->handle);
    return $file->handle;
}
#------------------------------------------------------------------------
=head2 close_handle

 Title   : close_handle
 Usage   : $obj->close_handle;
 Function: Closes open file handle in object. Should before creating
           threads or forking.
 Returns : True on success
 Args    : None

=cut
sub close_handle {
    my $self = $_[0];
    my $FH = $self->__file->handle;
    $self->__file->handle = undef; #clear value
    return close($FH) if($FH);
}
#------------------------------------------------------------------------
sub __open_handle {
    my $self = $_[0];
    my $file = $self->__file;

    my $arg = $file->abs_path;
    die "ERROR: File $arg does not exist" if(! -f $arg);
    sysopen(my $IN, $arg, O_RDONLY|O_BINARY) or die "ERROR: Could not open input file $arg: $!";
    sysseek($IN, $file->pos+length($file->io_buffer), SEEK_SET); #move to expected position
    $file->handle = $IN;

    return 1;
}
#------------------------------------------------------------------------
=head2 reset_handle

 Title   : close_handle
 Usage   : $obj->reset_handle;
 Function: Close and reopen file handle in object. Use after creating
           threads or forking if close_handle was not called.
 Returns : True on success
 Args    : None

=cut
sub reset_handle {
    return ($_[0]->close_handle && $_[0]->__open_handle) ? 1 : 0;
}
#------------------------------------------------------------------------
sub __io_buffer {
    return \ ($_[0]->__file->io_buffer); #always return ref
}
#------------------------------------------------------------------------
=head2 fill_io_buffer

 Title   : fill_io_buffer
 Usage   : $obj->fill_io_buffer($bytes);
 Function: Fill IO buffer for reading with given number of bytes. Methods
           like next_block read from this buffer. Usful to avoid multiple
           read operations if you expect to process all blocks in a large
           segment of a file.
 Returns : Size of IO buffer after fill
 Args    : Number of bytes (defaults to 65536 or 64kb when empty)

=cut
sub fill_io_buffer {
    my $self = $_[0];
    my $to_read = $_[1] || MAX_BSIZE;

    my $file = $self->__file;
    my $IN = $self->__handle;
    my $pos = $file->pos;
    my $size = $file->length;
    my $b_ref = $self->__io_buffer;

    return length($$b_ref) if($pos == $size);

    #read data into input buffer
    $pos = $pos + length($$b_ref); #position after accounting for buffer
    my $needed = ($size-$pos < $to_read) ? $size-$pos : $to_read;
    $self->__read($IN, $$b_ref, $needed, length($$b_ref))
	or die "ERROR: Problem reading file ".$file->path.": $!";

    return length($$b_ref); #return current buffer size
}
#------------------------------------------------------------------------
=head2 size_io_buffer

 Title   : size_io_buffer
 Usage   : $size = $obj->size_io_buffer();
 Function: Get current length of IO buffer
 Returns : Size of IO buffer
 Args    : None

=cut
sub size_io_buffer {
    return length(${$_[0]->__io_buffer});
}

#------------------------------------------------------------------------
#--------------      INFLATION BUFFER RELATED METHODS      --------------
#------------------------------------------------------------------------

sub __vbuffer_struct {
    return BGZFast::vbuffer->new(DATA => '',
				 UPSTREAM => undef,
				 LENGTH => 0,
				 POS => 0,
				 BLOCKS => [],
				 FEATURE_OK => undef);
}
#------------------------------------------------------------------------
sub __vbuffer {
    my $self = $_[0];
    $self->{VBUFFER} ||= $self->__vbuffer_struct();
    return $self->{VBUFFER};
}
#------------------------------------------------------------------------
=head2 vbuffer

 Title   : vbuffer
 Usage   : my $vbuffer = $obj->vbuffer()
 Function: Provides access to the virtual buffer data string.
 Returns : A string reference.
 Args    : None

=cut
sub vbuffer {
    my $self = $_[0];
    $self->__inflate_to_vbuffer(); #add block to buffer if present (also fixes buffer)
    return \ ($self->__vbuffer->data); #always return ref
}
#------------------------------------------------------------------------
=head2 vbuffer_offset

 Title   : vbuffer_offset
 Usage   : my $offset = $obj->vbuffer_offset(); #get current offset
           $obj->vbuffer_offset($offset); #move to given offset
           $obj->vbuffer_offset(0, SEEK_SET); #start of buffer
           $obj->vbuffer_offset(0, SEEK_END); #end of buffer
           $obj->vbuffer_offset(0, SEEK_CUR); #keeps current position
 Function: Get/set the offset in the virtual buffer. Will automaticaly
           trim blocks off of the virtual buffer if offset is set more
           than 64kb from the start. See Fcntl for wence values.
 Returns : Final buffer offset
 Args    : Number (negative numbers count backwards just as with seek)
           Wence (optional. values: SEEK_SET, SEEK_CUR, SEEK_END)

=cut

sub vbuffer_offset {
    my $self = $_[0];
    return $self->__vbuffer->pos if(!defined($_[1]));

    my $pos = $_[1];
    my $wence = $_[2] || SEEK_SET;
    my $buffer = $self->__vbuffer;

    #clear values made invalid by seek
    $self->__reset_seek_sensitive();

    #adjust position for wence
    my $go;
    my $len = length($buffer->data);
    if($wence == SEEK_SET){
	$go = $pos;
    }
    elsif($wence == SEEK_CUR){
	$go = $buffer->pos + $pos;
    }
    elsif($wence == SEEK_END){
	$go = $len + $pos;
    }
    else{
	die "ERROR: Invalid WENCE for seek operation";
    }
    
    if($go <= 0){
	$go = "0 but true"; #same behavior as sysseek
    }
    elsif($go > $len){
	$go = $len;
    }
    $buffer->pos = $go;

    #trim distant blocks in vbuffer
    $buffer->feature_ok = undef; #cannot guarantee feature position
    $self->__trim_vbuffer;

    return $buffer->pos;
}
#------------------------------------------------------------------------
sub __trim_vbuffer {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;

    #trim distant blocks in vbuffer
    if($buffer->pos >= MAX_ISIZE){
	my $bblocks = $buffer->blocks;
	while(@$bblocks){
	    my $bl = $bblocks->[0];
	    if($buffer->pos - MAX_ISIZE >= $bl->buffer_offset + $bl->i_length){ #past current position by 64kb
		shift @$bblocks; #remove block
		if($bl->i_length){ #non-empty block
		    $_->buffer_offset -= $bl->i_length foreach(@$bblocks);
		    $buffer->upstream = substr(substr($buffer->data, 0, $bl->i_length, ''), -1, 1); #trims buffer
		    $buffer->pos -= $bl->i_length;
		    $buffer->length -= $bl->i_length;
		}
	    }
	    else{
		last;
	    }
	}
    }

    return;
}
#------------------------------------------------------------------------
sub __reset_vbuffer {
    my $self = $_[0];

    #alter values inside buffer so changes carry over to references
    my $buffer = $self->__vbuffer;
    my $default = $self->__vbuffer_struct();
    if(ref($buffer eq 'HASH')){
	foreach my $key (keys %$buffer){
	    if(ref($buffer->{$key}) eq 'ARRAY'){
		@{$buffer->{$key}} = @{$default->{$key}};
	    }
	    elsif(ref($buffer->{$key}) eq 'HASH'){
		%{$buffer->{$key}} = %{$default->{$key}};
	    }
	    else{
		$buffer->{$key} = $default->{$key};
	    }
	}
    }
    else{
	for(my $i =0; $i < @$buffer; $i++){
	    if(ref($buffer->[$i]) eq 'ARRAY'){
		@{$buffer->[$i]} = @{$default->[$i]};
	    }
	    elsif(ref($buffer->[$i]) eq 'HASH'){
		%{$buffer->[$i]} = %{$default->[$i]};
	    }
	    else{
		$buffer->[$i] = $default->[$i];
	    }
	}
    }
    
    $self->__reset_seek_sensitive; #always affected
}
#------------------------------------------------------------------------
sub __inflate_to_vbuffer {
    my $self = $_[0];
    my $trim = $_[1] || 0;

    my $buffer = $self->__vbuffer;

    #get waiting block
    my $block = $self->{BLOCK};
    $self->{BLOCK} = undef; #reset block
    return unless($block && $block->d_length); #bad block

    #make sure new block is contiguous
    my $bblocks = $buffer->blocks;
    if(@$bblocks){
	my $last = $bblocks->[-1];
	if($last->file_offset+$last->d_length != $block->file_offset){ #not contiguous
	    $self->__reset_vbuffer;
	}
    }
	
    #trim distant blocks in buffer
    $self->__trim_vbuffer if($trim);

    #inflate block to buffer
    $block->buffer_offset = $buffer->length;
    if($block->i_length != 0){ #not empty
	$self->__inflate(\ ($block->data), \ ($buffer->data));
	$buffer->length += $block->i_length;
    }
    push(@{$buffer->blocks}, $block);

    $self->__set_upstream if(!defined($buffer->upstream));

    return $block->i_length;
}
#------------------------------------------------------------------------
=head2 grow_vbuffer

 Title   : grow_vbuffer
 Usage   : $obj->grow_vbuffer();
 Function: Adds next block in file to the virtual buffer without altering
           current buffer or block offset.
 Returns : Number (length of virtual buffer)
 Args    : Number (Optional target length to add to virtual buffer)

=cut
sub grow_vbuffer {
    my $self = $_[0];
    my $count = $_[1] || 1;

    my $added = 0;
    while($added < $count){
	last unless($self->__read_block());
	$added += $self->__inflate_to_vbuffer();
    }

    return $added;
}
#------------------------------------------------------------------------
#--------------            INDEX RELATED METHODS           --------------
#------------------------------------------------------------------------
sub __gzi_struct {
    return BGZFast::gzi->new(PATH     => undef,
			     ABS_PATH => undef,
			     LENGTH   => 0,
			     INDICES  => {BLOCKS => []});
}
#------------------------------------------------------------------------
sub __gzi {
    my $self = $_[0];
    $self->{GZI} ||= $self->__gzi_struct();
    return $self->{GZI};
}
#------------------------------------------------------------------------
=head2 gzi

 Title   : gzi
 Usage   : $obj->gzi();
 Function: Get/build the GZI index object. Builds and parses index if object
           does not yet exist.
 Returns : GZI index object.
 Args    : None

=cut
sub idx {return gzi(@_)} #sugar
sub gzi {
    my $self = $_[0];
    if(!$self->{GZI}){
	$self->{GZI} = (-f $self->gzi_file) ? $self->parse_gzi() : $self->make_gzi();
    }
    return $self->{GZI};
}
#------------------------------------------------------------------------
=head2 gzi_file

 Title   : gzi_file
 Usage   : $path = $obj->gzi_file();
           $path = $obj->gzi_file($file);
 Function: Get expected/observed path to GZI index file
 Returns : Path of GZI index file
 Args    : Path (optional path of file to GZI index. defaults to current file)

=cut
sub index_file {return(gzi_file(@_))} #sugar
sub gzi_file {
    my $self = $_[0];
    my $file = $_[1];

    if($file){
	return $file.".gzi";
    }
    elsif($self->{GZI} && defined($self->{GZI}->path)){
        return $self->{GZI}->path;
    }
    return $self->file.".gzi";
}
#------------------------------------------------------------------------
=head2 make_gzi

 Title   : make_gzi
 Usage   : $path = $obj->make_gzi();
 Function: Make GZI index of file. Calls to external indexer. Can Parallelize
           if CPUS specified in parameters given to 'new'.
 Returns : Path of GZI index file (undef on failure)
 Args    : None

=cut
sub make_gzi {
    my $self = $_[0];
    my $cpus = $self->{CPUS} || 1;

    my $file = $self->file;
    my $gzi_file = $self->gzi_file($file);
    if(my $exe = File::Which::which('pbgzip') || File::Which::which('bgzip')){
	#build command
	my @cmd = ($exe, '-i', '-f');
	if($exe =~ /pbgzip$/){
	    push(@cmd, '-p' => $cpus);
	}
	elsif($exe =~ /bgzip$/){
	    push(@cmd, '-@' => $cpus);
	}
	push(@cmd, $file);
	
	#run indexer
	my $stat = system(@cmd);
	unlink($gzi_file) if($stat); #on failure
    }    
    die "ERROR: ".__PACKAGE__." could not generate GZI index files" if(!-f $gzi_file);

    return $self->parse_gzi($gzi_file);
}
#------------------------------------------------------------------------
=head2 parse_gzi

 Title   : parse_gzi
 Usage   : $obj->parse_gzi();
           $obj->parse_gzi($gzi_file);
 Function: Parse given GZI index of file.
 Returns : True on success
 Args    : Path (optional path of GZI index. defaults to $obj->gzi_file)

=cut
sub parse_index {return parse_gzi(@_)}
sub parse_gzi {
    my $self = $_[0];
    my $gzi_file = $_[1] || $self->gzi_file;

    #not compressed
    my $size = (stat($gzi_file))[7];
    if(!$self->{IS_BGZF}){
	#hash representation
	my $gzi = $self->__gzi;
	$gzi->path = $gzi_file;
	$gzi->abs_path = Cwd::abs_path($gzi_file);
	$gzi->length = $size;
	$gzi->indices('BLOCKS') = [[0, $size, 0, $size]];
	return 1;
    }

    #read in entire file
    my $raw = ''; #raw data
    sysopen(my $IN, $gzi_file, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $gzi_file: $!";
    $self->__read($IN, $raw, $size, 0) or die "ERROR: Could not read $gzi_file: $!";
    close($IN);

    #first 8 bytes contain block entry count
    my $soffset = 0;
    my $nblocks = unpack('q<', substr($raw, $soffset, 8));
    $soffset += 8;
    die "ERROR: GZI index is empty" if($nblocks == -1);
    die "ERROR: GZI index appears to be corrupt" if($size != 8+$nblocks*16);

    #add block entries from index
    my @index = [0, undef, 0, undef]; #first block is always implicit in the index
    for(my $i = 1; $i <= $nblocks; $i++){
        my ($boffset, $ioffset) = unpack('q<q<', substr($raw, $soffset, 16));
        $soffset += 16;
        $index[-1][1] = $boffset - $index[-1][0];
        $index[-1][3] = $ioffset - $index[-1][2];
        push(@index, [$boffset, undef, $ioffset, undef]);
    }

    #last block is always incomplete and must be read from bgzf file
    my $filepath = $self->file;
    sysopen($IN, $filepath, O_RDONLY|O_BINARY) or die "ERROR: Could not open input file $filepath: $!";
    sysseek($IN, $index[-1][0], SEEK_SET);

    #determine remaining index values direct from bgzf file
    my $data = '';
    my $fsize = $self->__file->length;
    my $needed = $fsize - $index[-1][0] - 28; #expected terminal block size
    $self->__read($IN, $data, $needed, 0) or die "ERROR: Failure to read $filepath: $!";

    #add info for terminal blocks
    if(length($data)){
	#complete last block entry
	($index[-1][1],
	 $index[-1][3]) = $self->block_stat($data);

	#add an explicit EOF block index entry
	my $boffset = $index[-1][0]+$index[-1][1];
	my $ioffset = $index[-1][2]+$index[-1][3];
	push(@index, [$boffset, 28, $ioffset, 0]);
    }
    else{
	#complete last block entry
	$index[-1][1] = 28;
	$index[-1][3] = 0;
    }

    #verify final size matches index expectation
    if($index[-1][0]+$index[-1][1] != $fsize){
	die "ERROR: File size does not match index. GZI file may be corrupt";
    }
    
    #hash representation
    my $gzi = $self->__gzi;
    $gzi->path = $gzi_file;
    $gzi->abs_path = Cwd::abs_path($gzi_file);
    $gzi->length = $size;
    $gzi->indices('BLOCKS') = \@index;

    $self->{GZI} = $gzi;
    
    return $gzi;
}
#------------------------------------------------------------------------
sub __gzi2string {
    my $self = $_[0];
    my $gzi = $self->gzi;

    my $index = $gzi->indices('BLOCKS');

    #block count (ignore explicit first and last block entries)
    my $nblocks = scalar(@$index)-2 || -1;
    my $data = pack('q<', $nblocks);

    #add each index entry
    for(my $i = 1; $i < $#$index; $i++){
        $data .= pack('q<q<', $index->[$i][0], $index->[$i][2]);
    }

    return $data;
}
#------------------------------------------------------------------------
=head2 write_gzi

 Title   : write_gzi
 Usage   : $obj->write_gzi();
           $obj->write_gzi($file);
 Function: Write GZI index to a file.
 Returns : True on success
 Args    : Path (optional write path. defaults to $obj->gzi_file)

=cut
sub write_gzi {
    my $self = $_[0];
    my $gzi_file = $_[1] || $self->gzi_file();

    #write to file
    my $gzi_str = $self->__gzi2string();
    sysopen(my $OUT, $gzi_file, O_RDWR|O_CREAT|O_TRUNC)
        or die "ERROR: Could not open $gzi_file for writing: $!";
    $self->__write($OUT, $gzi_str) or die "ERROR: Failed writing to $gzi_file: $!";
    close($OUT);

    #update index file info
    my $gzi = $self->__gzi();
    $gzi->path = $gzi_file;
    $gzi->abs_path = Cwd::abs_path($gzi_file);
    $gzi->length = (stat($gzi_file))[7];

    return 1;
}
#------------------------------------------------------------------------
=head2 query_gzi

 Title   : query_gzi
 Usage   : ($real, $virtual) = $obj->query_gzi($pos);
 Function: Get real and virtual file offsets corresponding to the
           position query
 Returns : Real and virtual file offsets (just real in scalar context)
 Args    : Position to search GZI index for (virtual)

=cut
sub query_gzi {
    my $self = $_[0];
    my $pos = $_[1];

    return $pos if(!$self->{IS_BGZF});

    my $index = $self->gzi->indices('BLOCKS');
    foreach my $o (@$index){
	next unless($o->[2] <= $pos && $pos < $o->[2]+$o->[3]);
	return wantarray ? ($o->[0], $pos-$o->[0]) : $o->[0];
    }

    return wantarray ? ($self->__file->length, 0) : $self->__file->length;
}
#------------------------------------------------------------------------
#--------------            BLOCK RELATED METHODS           --------------
#------------------------------------------------------------------------

sub __block_struct {
    return BGZFast::block->new(DATA          => '',
			       I_LENGTH      => 0,
			       D_LENGTH      => 0,
			       BUFFER_OFFSET => undef,
			       FILE_OFFSET   => undef);
}
#------------------------------------------------------------------------
sub __block {
    my $self = $_[0];
    return $self->{BLOCK};
}
#------------------------------------------------------------------------
#scans inflation buffer for next block (assumes blocks and buffer already reset)
sub __adjust_to_next_block {
    my $self = $_[0];

    my $file = $self->__file;
    my $pos  = $file->pos;
    return $pos if($file->block_ok);

    #shortcut when not compressed
    if($pos == 0 || !$self->{IS_BGZF}) { #start or not compressed
	$file->block_ok = 1; #set flag
	return $pos;
    }

    #get first encountered block offset
    my $size = $file->length; #file size
    my $b_ref = $self->__io_buffer; #input buffer
    if($size-$pos < 28){ #too small for a block
	$pos = $self->__file_seek($size) if($pos != $size);
	$file->block_ok = 1;
	return $pos;
    }

    $self->fill_io_buffer(MAX_BSIZE) if(length($$b_ref) < MAX_BSIZE);
    my $offset = $self->find_block_offset($b_ref, 0);
    die "ERROR: Data is corrupt. Could not locate block header." if($offset == -1);
    $file->pos = $pos = $pos+$offset; #adjust file position with first block offset
    substr($$b_ref, 0, $offset, ''); #trim buffer for offset
    $file->block_ok = 1;

    return $pos;
}
#------------------------------------------------------------------------
#gets last character immediately upstream of current block
sub __set_upstream {
    my $self = $_[0];
    return if(!@{$self->__vbuffer->blocks}); #nothing to be upstream of

    my $IN = $self->__handle; #filehandle
    my $buffer = $self->__vbuffer;
    my $pos  = $buffer->blocks(0)->file_offset;
    my $file = $self->__file;
    my $b_ref = $self->__io_buffer;

    if($pos == 0) { #start of file
	$buffer->upstream = '';
	return;
    }
    elsif(!$self->{IS_BGZF}){ #not compressed
	#read upstream char from file
	sysseek($IN, $pos-1, SEEK_SET);
	$self->__read($IN, $buffer->upstream, 1) or
	    die "ERROR: Problem reading file ".$file->path.": $!";
	sysseek($IN, $file->pos+length($$b_ref), SEEK_SET); #restore expected
	return;
    }

    #get upstream char off of previous block
    my $up = '';
    my $raw = '';
    my $cur = $pos;
    my $alt = ($cur - MAX_BSIZE > 0) ? $cur - MAX_BSIZE : 0;
    while(1){
	#read in data
	my $to_read = $cur - $alt;
	my $copy = $raw;
	$raw = '';
	sysseek($IN, $alt, SEEK_SET);
	$self->__read($IN, $raw, $to_read, 0) or die "ERROR: Could not read from $file: $!";
	$raw .= $copy; #concatenate

	#reverse scan for block header
	my $offset = $self->rfind_block_offset(\$raw, length($raw)); #reverse scan
	if($offset == -1){
	    die "ERROR: Data is corrupt. Could not locate block header." if($alt == 0);
	    $cur = $alt;
	    $alt = ($cur - MAX_BSIZE > 0) ? $cur - MAX_BSIZE : 0;
	    next;
	}

	#get last character
	$self->__inflate(\ (substr($raw, $offset, length($raw)-$offset)), \$up); #replace required
	$up = substr($up, -1, 1); #last char
	last if(length($up) == 1 || ($alt == 0 && $offset == 0)); #loops on empty block
    }
    sysseek($IN, $file->pos+length($$b_ref), SEEK_SET); #restore expected file position
    $buffer->upstream = $up;
    return;
}
#------------------------------------------------------------------------
#assumes file position is already set to start of block
sub __read_block {
    my $self = $_[0];
    my $IN = $self->__handle;
    my $file = $self->__file;
    my $i_ref = $self->__io_buffer;

    #read data into input buffer
    $self->fill_io_buffer(MAX_BSIZE) if(length($$i_ref) < MAX_BSIZE);
    my $block = $self->{BLOCK} = $self->__block_from_ref($i_ref); #make block
    return 0 if(!$block);

    #fix values
    $block->file_offset += $file->pos;
    $file->pos += $block->d_length; #move file position
    substr($$i_ref, 0, $block->d_length, ''); #trim input buffer

    return 1;
}
#------------------------------------------------------------------------
#assumes file position is set to start of block
sub __block_from_ref {
    my $self = $_[0];
    my $ref = ref($_[1]) ? $_[1] : \$_[1];
    my $offset = $_[2] || 0;

    #not compreseed
    if(!$self->{IS_BGZF}){
	return undef if($offset == length($$ref)); #end of string
	my $block = $self->__block_struct();
        $block->data = substr($$ref, $offset, MAX_ISIZE); #uncompressed max is MAX_ISIZE
        $block->file_offset = $offset;
        $block->d_length = length($block->data);
        $block->i_length = length($block->data);
	return $block;
    }

    #read block from reference
    return undef if(length($$ref)-$offset < 28); #too short to contain block
    my $block = $self->__block_struct();
    ($block->d_length,
     $block->i_length) = $self->block_stat($ref, $offset);
    return undef if(!defined($block->d_length));
    $block->file_offset = $offset;
    $block->data = substr($$ref, $offset, $block->d_length);

    return $block;
}
#------------------------------------------------------------------------
=head2 current_block

 Title   : current_block
 Usage   : $block = $obj->current_block();
 Function: Get block overlapping current virtual buffer / block offsets
 Returns : block object
 Args    : None

=cut
sub current_block {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;

    $self->__inflate_to_vbuffer;
    return if(!@{$buffer->blocks});

    #find overlapping block
    my $block;
    my $pos = $buffer->pos;
    foreach my $bl (@{$buffer->blocks}){
	if($bl->buffer_offset <= $pos && $pos < $bl->buffer_offset + $bl->i_length){
	    $block = $bl;
	    last;
	}
    }

    #EOF block special case
    if(!$block){
	my $l = $buffer->blocks(-1);
	($block) = $l if($pos == $l->buffer_offset && $l->i_length == 0 &&
			 $l->file_offset + $l->d_length == $self->__file->length);
    }

    return $block;
}
#------------------------------------------------------------------------
=head2 next_block

 Title   : next_block
 Usage   : $block = $obj->next_block();
 Function: Get next block in file. Will move virtual buffer / block
           offsets to the start of this new block.
 Returns : block object
 Args    : None

=cut
sub next_block {
    my $self = $_[0];
    
    $self->__read_block();

    #fix buffer offset because of new block
    if(@{$self->__vbuffer->blocks} && $self->{BLOCK}){
	my $block = $self->__block();
	my $last = $self->__vbuffer->blocks(-1);
        if($last->file_offset+$last->d_length != $block->file_offset){ #not contiguous
            $self->__reset_vbuffer;
        }
	else{ #end of buffer will be start of this block when decompressed
	    $self->vbuffer_offset(0, SEEK_END);
	}
    }

    return $self->__block();
}
#------------------------------------------------------------------------
=head2 block_head

 Title   : block_head
 Usage   : $block = $obj->block_head();
 Function: Shatter the block overlapping the current virtual buffer /
           block offsets and return a new block of data upstream of the
           offsets.
 Returns : block object
 Args    : None

=cut
#temp
sub block_head {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;

    $self->__inflate_to_vbuffer(); #just in case
    my ($block) = $self->current_block();
    return '' if(!$block);

    my $offset = $block->buffer_offset;
    my $length = $buffer->pos - $block->buffer_offset;

    my $head = '';
    my $string = substr($buffer->data, $offset, $length);
    $self->__deflate(\$string, \$head) if(length($string));

    return $head;
}
#------------------------------------------------------------------------
=head2 block_tail

 Title   : block_head
 Usage   : $block = $obj->block_tail();
 Function: Shatter the block overlapping the current virtual buffer /
           block offsets and return a new block of data downstream of the
           offsets.
 Returns : block object
 Args    : None

=cut
sub block_tail {
    my $self = $_[0];
    my $buffer = $self->__vbuffer;

    $self->__inflate_to_vbuffer(); #just in case
    my ($block) = $self->current_block();

    #add another block if needed
    if(!$block){
	$self->__read_block();
	$self->__inflate_to_vbuffer();
	($block) = $self->current_block();
	return '' if(!$block); #end of file
    }

    my $offset = $buffer->pos;
    my $length = $block->i_length - ($buffer->pos - $block->buffer_offset);

    my $tail = '';
    my $string = substr($buffer->data, $offset, $length);
    $self->__deflate(\$string, \$tail) if(length($string));

    return $tail;
}
#------------------------------------------------------------------------
=head2 write_block

 Title   : write_block
 Usage   : $obj->write_block($FH);
 Function: Write a given block to a given file handle.
 Returns : Bytes written (undef on failure with $! set)
 Args    : File handle to write to
           Block object to write

=cut
sub write_block {
    my $self = $_[0];
    my $FH   = $_[1];
    my $block = $_[2];

    $self->__write($FH, $block->data) or die "ERROR: Failed writing to file: $!";
}
#------------------------------------------------------------------------
=head2 write_as_blocks

 Title   : write_as_blocks
 Usage   : $obj->write_as_blocks($FH);
 Function: Compress a given string into BGZF blocks and write them to a
           given file handle.
 Returns : Bytes written (undef on failure with $! set)
 Args    : File handle to write to
           String/String ref

=cut
sub write_as_blocks {
    my $self = $_[0];
    my $FH   = $_[1];
    my $ref = (ref($_[2])) ? $_[2] : \ ($_[2]); #make into reference
    
    $self->__write($FH, $self->zip($ref)) or die "ERROR: Failed writing to file: $!";
}
#------------------------------------------------------------------------
=head2 write_eof_block

 Title   : write_eof_block
 Usage   : $obj->write_eof_block($FH);
 Function: Write an EOF block to a given file handle.
 Returns : Bytes written (undef on failure with $! set)
 Args    : File handle to write to

=cut
sub write_eof_block {
    my $self = $_[0];
    my $FH   = $_[1];
    my $data = $self->eof_block;

    $self->__write($FH, $data) or die "ERROR: Failed writing to file: $!";
}
#------------------------------------------------------------------------
=head2 write_section

 Title   : write_section
 Usage   : $obj->write_section($file_path, $vstart, $vend);
 Function: Will write all data between two virtual position boundaries
           to the given file handle. Will trim blocks around start and
           end positions appropriately and add eof block. 
 Returns : Bytes written (undef on failure with $! set)
 Args    : File path to write to
           Start position
           End position

=cut
sub write_section {
    my $self  = $_[0];
    my $path  = $_[1];
    my $start = $_[2];
    my $end   = $_[3];

    #identify parts to copy
    $self = $self->clone(); #use copy to avoid buffer clearing
    my ($off_end, $vend) = $self->pos_seek($end);
    my $tail = $self->block_head;
    my ($off_beg, $vbeg) = $self->pos_seek($start);
    my $head = '';
    if($off_beg == $off_end){ #if starts and ends on same block
	#rebuild $tail
	my $string = '';
	$self->__inflate(\$tail, \$string) if(length($tail));
	$string = substr($string, -($vend-$vbeg));
	$self->__deflate(\$string, \$tail) if(length($string));
    }
    elsif($vbeg != 0){ #only if block is split
	$head .= $self->block_tail;
	$off_beg += $self->current_block->d_length; #shift passed first partial block
    }
    $tail .= $self->eof_block() if($self->{IS_BGZF});

    #copy data to output file
    sysopen(my $OUT, $path, O_RDWR|O_CREAT|O_TRUNC|O_BINARY)
	or die "ERROR: Could not open output file $path: $!";
    if(length($head)){
	$self->__write($OUT, $head) or die "ERROR: Could not write to $path: $!";
    }
    if($off_end - $off_beg){
	my $IN = $self->__handle;
	sysseek($IN, $off_beg, 0);
	my $needed = ($off_end - $off_beg);
	$self->__copy_data($IN, $OUT, $needed)
	    or die "ERROR: Could not copy data from ".$self->file." to $path: $!";
    }
    if(length($tail)){
	$self->__write($OUT, $tail) or die "ERROR: Could not write to $path: $!";
    }
    close($OUT);

    return length($head)+($off_end-$off_beg)+length($tail);
}

#------------------------------------------------------------------------
#--------------            DATA SEEK/TELL METHODS          --------------
#------------------------------------------------------------------------

sub __file_seek {
    my $self = $_[0];
    my $pos = $_[1];
    my $wence = $_[2] || SEEK_SET;

    my $file = $self->__file;
    my $IN = $self->__handle;

    #handle wence
    my $new;
    my $old = $file->pos;
    my $len = $file->length;
    if($wence == SEEK_SET){
        $new = $pos;
    }
    elsif($wence == SEEK_CUR){
        $new = $old + $pos;
    }
    elsif($wence == SEEK_END){
        $new = $len + $pos;
    }
    else{
        die "ERROR: Invalid WENCE for seek operation";
    }

    #correct pos
    if($new <= 0){
        $new = "0 but true"; #same behavior as sysseek
    }
    elsif($new > $len){
        $new = $len;
    }
    $file->pos = $new;

    #clear values made invalid by seek
    $self->{BLOCK} = undef;
    $self->__reset_vbuffer;
    $file->block_ok = undef; #unset flag that bock adjustment is set

    #adjust io read buffer
    my $b_ref = $self->__io_buffer;
    if($old < $new && $new < $old+length($$b_ref)){
	substr($$b_ref, 0, $new-$old, '');
	sysseek($IN, $new+length($$b_ref), SEEK_SET); #just after buffer
    }
    elsif($old != $new){
	$$b_ref = '';
	sysseek($IN, $new, SEEK_SET);
    }

    return $new;
}
#------------------------------------------------------------------------
sub __file_tell {
    return $_[0]->__file->pos;
}
#------------------------------------------------------------------------
=head2 block_seek

 Title   : block_seek
 Usage   : $obj->block_seek($reak, $virtual);
 Function: Move to a given real and virtual offset in a BGZF file.
 Returns : Real and virtual offset moved to (just real offset in scalar
           context)
 Args    : Real file offset
           Virtual file offset

=cut
sub block_seek {
    my $self = $_[0];
    my $file_pos = $_[1];
    my $block_pos = $_[2] || 0;

    my $buffer = $self->__vbuffer;

    #clear values made invalid by seek
    $self->__reset_seek_sensitive();

    #see if block is already active
    my $block;
    if(!$self->{IS_BGZF}){ #make positions absolute for uncompressed data
	my $abs_pos = $file_pos + $block_pos;
	($block) = grep {$_ && $_->file_offset <= $abs_pos &&
			 $abs_pos < $_->file_offset+$_->i_length} (@{$buffer->blocks},
								       $self->__block());
    }
    else{
	($block) = grep {$_ && $_->file_offset == $file_pos} (@{$buffer->blocks},
								$self->__block());
    }

    #go to new block if needed
    if(!$block){
	$self->__file_seek($file_pos); #move to before desired pos
	$self->__adjust_to_next_block(); #find block and set upstream char
	$self->__read_block();
	$block = $self->__block();
    }
    
    #set correct position in vbuffer
    die "ERROR: Position is outside of block" if($block_pos > $block->i_length);
    if(defined($block->buffer_offset)){
	$self->vbuffer_offset($block->buffer_offset + $block_pos, SEEK_SET);
    }
    elsif(@{$buffer->blocks} && $block){
        #make sure new block is contiguous
        my $last = $buffer->blocks(-1);
        if($last->file_offset+$last->d_length != $block->file_offset){ #not contiguous
            $self->__reset_vbuffer;
        }
	else{
	    $self->vbuffer_offset(0, SEEK_END);
	}
    }

    return $self->block_tell;
}
#------------------------------------------------------------------------
=head2 block_tell

 Title   : block_tell
 Usage   : ($real, $virtual) = $obj->block_seek()
 Function: Get current real and virtual offsets in a BGZF file.
 Returns : Real and virtual offsets (just real offset in scalar context)
 Args    : None

=cut
sub block_tell {
    my $self = $_[0];

    my $buffer = $self->__vbuffer;

    #need to be at a block or in a block
    $self->__adjust_to_next_block() if(!$self->__file->block_ok);

    #get needed values from corresponding block
    my ($block) = $self->current_block();
    if(!$block){
	return wantarray ? ($self->__file_tell, 0) : $self->__file_tell;
    }

    #return absolute offset for uncompressed data
    if(!$self->{IS_BGZF}){
	my $abs_offset = $block->file_offset + $buffer->pos - $block->buffer_offset;
	return wantarray ? ($abs_offset, 0) : $abs_offset;
    }

    #return file and block offsets
    return wantarray ?
	($block->file_offset, $buffer->pos - $block->buffer_offset) : $block->file_offset;
}
#------------------------------------------------------------------------
=head2 smart_seek

 Title   : block_tell
 Usage   : ($real, $virtual) = $obj->smart_seek($pos)
           ($real, $virtual) = $obj->smart_seek($pos, SEEK_SET);
 Function: Move to a given position and auto adjust to start of the next
           BGZF block in the file. See Fcntl for wence options.
 Returns : Real and virtual offsets (just real offset in scalar context)
 Args    : Number (negative will count backwards just lie seek)
           Wence (optional. values: SEEK_SET, SEEK_CUR, SEEK_END)

=cut
sub smart_seek {
    my $self = $_[0];
    my $pos = $_[1];
    my $wence = $_[2] || SEEK_SET;

    my $file = $self->__file;

    #handle wence
    my $new;
    my $old = $file->pos;
    my $len = $file->length;
    if($wence == SEEK_SET){
        $new = $pos;
    }
    elsif($wence == SEEK_CUR){
        $new = $old + $pos;
    }
    elsif($wence == SEEK_END){
        $new = $len + $pos;
    }
    else{
        die "ERROR: Invalid WENCE for seek operation";
    }
    $self->__file_seek($new); #move to before desired pos
    $self->__adjust_to_next_block();

    return $self->block_tell();
}

#------------------------------------------------------------------------
#--------------            BGZF_UTILITY ALIASES            --------------
#------------------------------------------------------------------------

#inflate/deflate compression block
sub __inflate { shift @_; return BGZFast::bgzf_utility::__inflate(@_); }
sub __deflate { shift @_; return BGZFast::bgzf_utility::__deflate(@_); }

#find block offsets in a scalar ref
sub validate_block_header {shift @_; return BGZFast::bgzf_utility::validate_block_header(@_); }
sub block_stat {shift @_; return BGZFast::bgzf_utility::block_stat(@_); }
sub find_block_offset {shift @_; return BGZFast::bgzf_utility::find_block_offset(@_); }
sub find_full_block {shift @_; return BGZFast::bgzf_utility::find_full_block(@_); }
sub rfind_block_offset {shift @_; return BGZFast::bgzf_utility::rfind_block_offset(@_); }
sub rfind_full_block {shift @_; return BGZFast::bgzf_utility::rfind_full_block(@_); }

#bgzip all blocks from input scalar ref
sub zip { shift @_;return BGZFast::bgzf_utility::zip(@_); }

#unbgzip all blocks from input scalar ref
sub unzip { shift @_; return BGZFast::bgzf_utility::unzip(@_); }
sub unzip_full { return BGZFast::bgzf_utility::unzip($_[1], $_[2], 1); }
sub unzip_partial { return BGZFast::bgzf_utility::unzip($_[1], $_[2], 0); }

#convenience method
sub sparse_check { shift @_;return BGZFast::bgzf_utility::sparse_check(@_); }

#end of file signature block
sub max_bsize { return MAX_BSIZE; }
sub max_isize { return MAX_ISIZE; }
sub eof_block { return EOF_BLOCK; }

#------------------------------------------------------------------------
#--------------      HOUSEKEEPING/CONVENIENCE METHODS      --------------
#------------------------------------------------------------------------
sub __write {
    my $self   = $_[0];
    my $OUT    = $_[1];
    #data is in $_[2]
    my $length = $_[3] || length($_[2]);
    my $offset = $_[4] || 0;

    my $needed = $length;
    while($needed){
	my $to_write = (33554432 > $needed) ? 33554432 : $needed;
	my $stat = syswrite($OUT, $_[2], $to_write, $offset+$length-$needed);
	return undef if(!defined($stat));
	$needed -= $stat;
    }

    return $length-$needed || "0 but true";
}

sub __read {
    my $self   = $_[0];
    my $IN     = $_[1];
    #data is in $_[2]
    my $length = $_[3] || (stat($IN))[7]-sysseek($IN, 0, SEEK_CUR);
    my $offset = $_[4] || 0;

    my $needed = $length;
    while($needed){
	my $to_read = ($needed > 33554432) ? 33554432 : $needed;
	my $stat = sysread($IN, $_[2], $to_read, $offset+$length-$needed);
	return undef if(!defined($stat));
	last if($stat == 0); #end of file
	$needed -= $stat;
    }

    return $length-$needed || "0 but true";
}

sub __copy_data {
    my $self = $_[0];
    my $IN   = $_[1];
    my $OUT  = $_[2];
    my $length = $_[3] || (stat($IN))[7]-sysseek($IN, 0, SEEK_CUR);

    my $rneeded = $length;
    while($rneeded){
	#read from infile
	my $data = '';
        my $to_read = ($rneeded > 33554432) ? 33554432 : $rneeded;
        my $rstat = sysread($IN, $data, $rneeded, 0);
        return undef if(!defined($rstat));
	last if($rstat == 0); #end of file
        $rneeded -= $rstat;

	#write same amount read to outfile
	my $wneeded = $rstat;
	while($wneeded){
	    my $wstat = syswrite($OUT, $data, $rstat, $rstat-$wneeded);
	    return undef if(!defined($wstat));
	    $wneeded -= $wstat;
	}
    }

    return $length-$rneeded || "0 but true";
}

#------------------------------------------------------------------------
#--------------                 FUNCTIONS                  --------------
#------------------------------------------------------------------------

sub _cap_hash {
    my $rhash = $_[0];
    my %hash = map {
        my $k = $_;
	my $v = $rhash->{$k};
	$k =~ s/^\-{1,2}//;
        $k =~ tr/a-z/A-Z/;
	$k => $v;
    } keys(%{$rhash});
    return \%hash;
}

1;
