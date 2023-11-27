#------------------------------------------------------------------------
#----                     UGP::vcf_manipulator                       ----
#------------------------------------------------------------------------
package UGP::vcf_manipulator;

use strict;
use warnings;
use Compress::Zlib;
use Fcntl qw(:DEFAULT :seek); #SEEK_SET, SEEK_CUR, SEEK_END

our $VERSION = '0.01';

#------------------------------------------------------------------------
#--------------       OBJECT INITIALIZATION METHODS        --------------
#------------------------------------------------------------------------

sub new {
    my $class = (ref($_[0])) ? ref(shift @_) : shift @_; #No copy constructor
    my $self = {};
    bless($self, $class);

    #initialize
    $self->__set_index_param;
    $self->init(@_);

    return $self;
}

sub init {
    my $self = shift;
    my @args = @_;

    #get user supplied parameters
    my $param;
    if (ref($args[0]) eq 'HASH') {
        $param = _cap_hash($args[0]);
    } else {
	@args = (FILE => $args[0]) if(@args == 1);
        $param = _cap_hash({@args});
    }

    #compression setting
    my $bgzf;
    if(defined($param->{BGZIP})){
	 $self->{BGZIP} = $bgzf = delete($param->{BGZIP});
    }

    #get file parameter
    if(defined(my $arg = delete($param->{FILE}))){
	die "ERROR: File $arg does not exist\n" if(! -f $arg);

	#detect compression
	$self->{BGZIP} = 1 if($arg =~ /\.gz$/ && !defined($bgzf)); #BGZF extension

	#get file
	sysopen(my $IN, $arg, O_RDONLY|O_BINARY) or die "ERROR: Could not open input file $arg: $!";
	if($self->{BGZIP}){
	    sysseek($IN, -28, SEEK_END); #last 28 bytes
	    my $last;
	    my $stat = sysread($IN, $last, 28, 0);
	    die "ERROR: The file $arg is missing the EOF marker\n"
		if($stat != 28 || $last ne $self->eof_block());
	    sysseek($IN, 0, SEEK_SET); #go back to start of file
	}

	#set file info
	my $file = $self->__file;
	$file->{PATH}     = $arg;
	($file->{NAME})   = $arg =~ /([^\/]+)$/;
	$file->{ABS_PATH} = Cwd::abs_path($arg);
	$file->{HANDLE}   = $IN;
	$file->{POS}      = 0;
	$file->{LENGTH}   = (stat($IN))[7];
    }
    else{
	die "ERROR: You must provide a file to initialize ".__PACKAGE__."\n";
    }

    #check for unrecognized parameters
    die "ERROR: Invalid option(s) (".join(', ', keys %$param).") passed to ".__PACKAGE__."::init\n"
	if(keys %$param);

    $self->__parse_header();

    return 1;
}

#------------------------------------------------------------------------
#--------------    OBJECT SPECIFIC METHODS TO OVERRIDE     --------------
#------------------------------------------------------------------------

#sets index parameters for specific file type (should be overridden by child module)
sub __set_index_param {
    my $self = shift;

    $self->{IDX_TYPE} = 2; #0: generic (1 based); 1: SAM; 2: VCF, 65536: generic (zero based)
    $self->{IDX_COL} = [0, 1, 1]; #column for seq, beg, end
    $self->{IDX_META} = '#'; #comment character
    $self->{IDX_SKIP} = 0; #number of lines to skip at beginning
    
    return;
}

#parses info from feature text (should be overridden by child module)
sub __fill_feature_ref {
    my $self = shift;
    my $feat = shift;

    my @F = split(/\t/, $feat->{DATA});
    chomp($F[-1]);

    $feat->{REFID} = $self->name2seqid($F[0]);
    $feat->{POS} = $F[1]-1; #make zero based
    $feat->{END} = $F[1]+length($F[3])-1; #make zero based
    $feat->{NEXT_REFID} = -1;
    $feat->{NEXT_POS} = -1;
    
    #add breakpoints for SVs
    if($F[4] =~ /[\[\]]\s*([^\[\]\:\s]+)\s*\:\s*([^\[\]\:\s]+)\s*[\[\]]/){
	$feat->{NEXT_POS} = $2-1; #make zero based
	$feat->{NEXT_REFID} = $self->name2seqid($1);
    }

    return;
}

#parses header reference text (should be overridden by child module)
sub __fill_header_ref {
    my $self = shift;
    my $header = shift;
    my $text = $header->{TEXT};

    #parse reference lines
    my $i = 0;
    foreach my $line (split(/\n/, $text)){
	next unless($line =~ /^##contig=/);
	my ($name) = $line =~ /[\,\<]ID=([^\,\>]+)/;
	my ($l_ref) = $line =~ /[\,\<]length=([^\,\>]+)/;
	
	$header->{REF}[$i]{NAME} = $name;
	$header->{REF}[$i]{LENGTH} = $l_ref;
	$header->{NAME2SEQID}{$name} = $i;
	$header->{NAME2LENGTH}{$name} = $l_ref;
	$i++;
    }
    $header->{N_REF} = $i;

    return;
}

#------------------------------------------------------------------------
#--------------              GENERIC METHODS               --------------
#------------------------------------------------------------------------

sub __file {
    my $self = shift;
    $self->{FILE} ||= {NAME     => undef,
		       PATH     => undef,
		       ABS_PATH => undef,
		       HANDLE   => undef,
		       LENGTH   => 0,
		       POS      => 0};
    return $self->{FILE};
}

sub __handle {
    my $self = shift;
    return $self->__file->{HANDLE};
}

sub __block {
    my $self = shift;
    $self->{BLOCK} ||= {DATA          => '',
			I_LENGTH      => 0,
			D_LENGTH      => 0,
			BUFFER_OFFSET => undef,
			FILE_OFFSET   => undef};
    return $self->{BLOCK};
}

sub __buffer {
    my $self = shift;
    $self->{BUFFER} ||= {DATA => '',
			 UPSTREAM => '',
			 LENGTH => 0,
			 POS => 0,
			 BLOCKS => []};
    return $self->{BUFFER};
}

sub __reset_buffer {
    my $self = shift;
    my $buffer = $self->__buffer;

    %$buffer = (DATA => '',
		UPSTREAM => '',
		LENGTH => 0,
		POS => 0,
		BLOCKS => []);
}

sub __feature {
    my $self = shift;

    $self->{FEATURE} ||= {DATA => '',
			  LENGTH => 0,
			  FILE_OFFSET => undef,
			  BLOCK_OFFSET => undef,
			  REFID => undef,
			  POS => undef,
			  END => undef,
			  NEXT_REFID => undef, #only on SV
			  NEXT_POS => undef}; #only on SV

    return $self->{FEATURE};
}

sub __header {
    my $self = shift;

    $self->{HEADER} ||= {DATA => '',
			 TAIL => '',
			 TAIL_FILE_OFFSET => undef,
			 TAIL_BLOCK_OFFSET => undef,
			 D_LENGTH => 0,
			 TEXT => '',
			 N_REF => undef,
			 NAME2SEQID => {},
			 NAME2LENGTH => {},
			 REF => []}; #REF has NAME and LENGTH

    return $self->{HEADER};
}

sub __idx {
    my $self = shift;
    $self->{IDX} ||= {PATH     => undef,
		      ABS_PATH => undef,
		      TYPE => $self->{IDX_TYPE}, #0: generic (1 based); 1: SAM; 2: VCF, 65536: generic (zero based)
		      COL  => $self->{IDX_COL},  #column for seq, beg, end
		      META => $self->{IDX_META}, #comment character
		      SKIP => $self->{IDX_SKIP}, #number of lines to skip at beginning
		      LENGTH   => 0,
		      INDICES => [],
		      N_NO_COOR => undef,
		      NAMES => [],
		      NAME2TBID => {},
		      SEQID2TBID => [],
		      TBID2SEQID => [],
		      REF => []};

    return $self->{IDX};
}

sub header {
    my $self = shift;

    return $self->{HEADER} if($self->{HEADER});

    $self->__parse_header;
    return $self->{HEADER};
}

sub idx {
    my $self = shift;

    if(!$self->{IDX}){
	$self->make_index();
	$self->__parse_idx();
    }

    return $self->{IDX};
}

sub tbi { return shift->idx}
sub idx_file {return shift->file.".tbi"}
sub make_index {
    my $self = shift;
    my $file = shift || $self->file;
    my $force = shift;

    my $idx_file = $self->idx_file;
    unlink($idx_file) if($force);
    if(! -f $idx_file){
	if(my $exe = File::Which::which('tabix')){
	    #build command for index type => 0: generic (1 based); 1: SAM; 2: VCF, 65536: generic (zero based)
	    my @cmd = ($exe);
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
	    unlink($idx_file) if($stat); #on failure
	}
    }
    die "ERROR: ".__PACKAGE__." could not generate TBI index files\n" if(!-f $idx_file);
}

#get/set vcf file
sub file {
    return shift->__file->{PATH};
}

sub __close_handle {
    return close(delete(shift->__file->{HANDLE}));
}

sub __reopen_handle {
    my $self = shift;
    my $file = $self->__file;

    my $arg = $file->{ABS_PATH};
    die "ERROR: File $arg does not exist\n" if(! -f $arg);
    sysopen(my $IN, $arg, O_RDONLY|O_BINARY) or die "ERROR: Could not open input file $arg: $!";
    sysseek($IN, $file->{POS}, SEEK_SET); #move to proper position
    $file->{HANDLE}   = $IN;

    return 1;
}

sub file_seek {
    my $self = shift;
    my $pos = shift;
    my $wence = shift || SEEK_SET;
    my $file = $self->__file;
    my $IN = $self->__handle;

    my $stat = sysseek($IN, $pos, $wence);
    $file->{POS} = $stat;

    #clear values made invalid by seek
    $self->{FEATURE} = undef  if($self->{FEATURE});
    $self->{BLOCK} = undef if($self->{BLOCK});
    $self->__reset_buffer if($self->{BUFFER});

    return $stat;
}

sub file_tell {
    return shift->__file->{POS};
}

sub buffer_seek {
    my $self = shift;
    my $pos = shift;
    my $wence = shift || SEEK_SET;
    my $buffer = $self->__buffer;

    #clear values made invalid by seek
    $self->{FEATURE} = undef  if($self->{FEATURE});

    my $go;
    my $len = length($buffer->{DATA});
    if($wence == SEEK_SET){
	$go = $pos;
    }
    elsif($wence == SEEK_SET){
	$go = $buffer->{POS} + $pos;
    }
    elsif($wence == SEEK_END){
	$go = $len + $pos;
    }
    else{
	die "ERROR: Invalid WENCE for seek operation\n";
    }

    if($go <= 0){
	$go = "0 but true"; #same behavior as sysseek
    }
    elsif($go > $len){
	$go = $len;
    }
    $buffer->{POS} = $go;

    return $go;
}

sub buffer_tell {
    return shift->__buffer->{POS};
}

sub block_seek {
    my $self = shift;
    my $file_pos = shift;
    my $block_pos = shift;
    my $buffer = $self->__buffer;

    #clear values made invalid by seek
    $self->{FEATURE} = undef  if($self->{FEATURE});

    #see if block is already active
    my $block;
    $self->__inflate_to_buffer;
    if(!$self->{BGZIP}){ #make positions absolute for uncompressed data
	my $abs_pos = $file_pos + $block_pos;
	($block) = grep {$_->{FILE_OFFSET} <= $abs_pos &&
			     $abs_pos < $_->{FILE_OFFSET}+$_->{I_LENGTH}} @{$buffer->{BLOCKS}};
    }
    else{
	($block) = grep {$_->{FILE_OFFSET} == $file_pos} @{$buffer->{BLOCKS}};
    }

    #go to new block if needed
    if(!$block){
	$self->file_seek($file_pos);
	$self->__adjust_to_next_block();
	$self->__read_block();
	$self->__inflate_to_buffer();
	$block = $buffer->{BLOCKS}[0];
    }
    die "ERROR: Position is outside of block\n" if($block_pos > $block->{I_LENGTH});
    $buffer->{POS} = $block->{BUFFER_OFFSET} + $block_pos;

    return $self->block_tell;
}

sub block_tell {
    my $self = shift;
    my $buffer = $self->__buffer;

    #get needed values from corresponding block
    my $all = $buffer->{BLOCKS};
    my ($block) = grep {$buffer->{POS} < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$all};
    ($block) = grep {$buffer->{POS} == $_->{BUFFER_OFFSET} && $_->{I_LENGTH} == 0} @{$all} if(!$block); #EOF block

    #return absolute offset for uncompressed data
    if(!$self->{BGZIP}){
	($block) = grep {$buffer->{POS} == $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$all} if(!$block); #last block
	my $abs_offset = $block->{FILE_OFFSET} + $buffer->{POS} - $block->{BUFFER_OFFSET};
	return wantarray ? ($abs_offset, 0) : $abs_offset;
    }

    #return file and block offsets
    return wantarray ?
	($block->{FILE_OFFSET}, $buffer->{POS} - $block->{BUFFER_OFFSET}) : $block->{FILE_OFFSET};
}

#jump to a given byte position in the file (will auto-adjust to start of next read in block)
sub smart_seek {
    my $self = shift;
    my $file = $self->__file;
    my $header = $self->header;

    my $pos = $self->file_seek(@_);
    $self->file_seek($file->{LENGTH}-28) if($self->{BGZIP} && $pos >= $file->{LENGTH}-28); #jump back from eof block
    $self->file_seek($header->{TAIL_FILE_OFFSET}) if($pos < $header->{TAIL_FILE_OFFSET});
    $self->__adjust_to_next_block();

    #keep reading blocks until feature is found
    while($self->__read_block()){
	$self->__inflate_to_buffer();
	last if(defined($self->__adjust_to_next_feature()));
    }

    return $self->block_tell;
}

#seek to an exact position
sub pos_seek {
    my $self = shift;
    my $chr = shift;
    my $pos = shift;
    my $header = $self->header;
    my $tbi    = $self->tbi;

    my $offset;
    my $v_offset;
    my $seqid = $self->name2seqid($chr);
    if($seqid == -1){
	#return offset of eof_block if seeking end of unmapped bin
	if($pos >= 0){
	    return wantarray ? ($self->__file->{LENGTH}-28, 0) : $self->__file->{LENGTH}-28;
	}

	#set start to be end of last mapped bin
	for(my $i = $#{$tbi->{INDICES}}; $i >= 0; $i--){
	    my $index = $tbi->{INDICES}[$i]{LINEAR};
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
	for(my $i = $seqid; $i < $header->{N_REF}; $i++){
	    my $tbid = $tbi->{SEQID2TBID}[$i]; #convert to tabix ID
	    next if(!defined($tbid));

	    my $index = $tbi->{INDICES}[$tbid]{LINEAR};
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
	    for(my $i = $#{$tbi->{INDICES}}; $i >= 0; $i--){
		my $index = $tbi->{INDICES}[$i]{LINEAR};
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
		return wantarray ? ($self->__file->{LENGTH}-28, 0) : $self->__file->{LENGTH}-28;
	    }
	}
    }

    #when the index is empty (no mapped reads) use the end of header position
    if(!defined($offset) && !defined($v_offset)){
	($offset, $v_offset) = ($header->{TAIL_FILE_OFFSET}, $header->{TAIL_BLOCK_OFFSET});
    }

    #go to desired position and then adjust for feature match
    $self->block_seek($offset, $v_offset);
    while(my $feat = $self->next_feature){
	my $seqid2 = $feat->{REFID};
	my $pos2 = $feat->{POS}; #zero based
	next if($seqid2 < $seqid && $seqid2 != -1); #too far in front
	next if($seqid2 > $seqid && $seqid == -1); #too far in front
	next if($seqid2 == $seqid && $pos2 < $pos); #too far in front
	return $self->block_seek($feat->{FILE_OFFSET}, $feat->{BLOCK_OFFSET});
    }

    #return offset of eof_block
    return wantarray ? ($self->__file->{LENGTH}-28, 0) : $self->__file->{LENGTH}-28;
}

sub current_block {
    my $self = shift;
    my $buffer = $self->__buffer;

    $self->__inflate_to_buffer;
    my ($block) = grep {$buffer->{POS} < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$buffer->{BLOCKS}};
    ($block) = grep {$buffer->{POS} == $_->{BUFFER_OFFSET} && $_->{I_LENGTH} == 0} @{$buffer->{BLOCKS}} if(!$block); #EOF block

    return $block;
}

sub next_block {
    my $self = shift;

    $self->__read_block();

    return $self->__block;
}

#write from an existing block
sub write_block {
    my $self = shift;
    my $FH   = shift;
    my $block = shift || $self->__block;

    syswrite($FH, $block->{DATA});
}

#create block from string then write it
sub write_to_block {
    my $self = shift;
    my $FH   = shift;
    my $string = (ref($_[0])) ? shift @_ : \ ($_[0]); #make into reference
    
    my $data;
    $self->__deflate($string, \$data);
    syswrite($FH, $data);
}

sub write_eof_block {
    my $self = shift;
    my $FH   = shift;
    my $data = $self->eof_block;

    syswrite($FH, $data);
}

sub block_head {
    my $self = shift;
    my $buffer = $self->__buffer;

    my ($block) = grep {$buffer->{POS} < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$buffer->{BLOCKS}};

    my $offset = ($block) ? $block->{BUFFER_OFFSET} : 0;
    my $length = ($block) ? $buffer->{POS} - $block->{BUFFER_OFFSET} : 0;

    my $head = '';
    my $string = substr($buffer->{DATA}, $offset, $length);
    $self->__deflate(\$string, \$head) if(length($string));

    return $head;
}

sub block_tail {
    my $self = shift;
    my $buffer = $self->__buffer;

    my ($block) = grep {$buffer->{POS} < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$buffer->{BLOCKS}};

    #add another block if needed
    if(!$block){
	$self->__inflate_to_buffer(); #just in case on is already loaded
	$self->__read_block();
	$self->__inflate_to_buffer();
	($block) = grep {$buffer->{POS} < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @{$buffer->{BLOCKS}};
	return '' if(!$block); #end of file
    }

    my $offset = $buffer->{POS};
    my $length = ($block) ? $block->{I_LENGTH} - ($buffer->{POS} - $block->{BUFFER_OFFSET}) : $block->{I_LENGTH};

    my $tail = '';
    my $string = substr($buffer->{DATA}, $offset, $length);
    $self->__deflate(\$string, \$tail) if(length($string));

    return $tail;
}

sub next_feature {
    my $self = shift;
    my $buffer = $self->__buffer;
    my $feat = $self->__feature;

    #adjust to position if not yet set properly
    my $upstream = ($buffer->{POS} == 0) ?
        $buffer->{UPSTREAM} : substr($buffer->{DATA}, $buffer->{POS}-1, 1);
    if($buffer->{UPSTREAM} ne "\n"){
	$self->__adjust_to_next_feature();
    }

    #process file feature
    my $meta = $self->{IDX_META};
    my $data_len = $buffer->{LENGTH};
    my $block_offset = $buffer->{POS};
    while(1){
	my $stat;
	while(!($stat = index($buffer->{DATA}, "\n", $block_offset)+1)){
	    return unless($self->__read_block);
	    my $block = $self->__block;
	    my $e_len = $buffer->{LENGTH} + $block->{I_LENGTH}; #expected length
	    $self->__inflate_to_buffer;
	    $data_len = $buffer->{LENGTH};
	    $block_offset -= $e_len - $data_len; #accounts for potential truncation
	}
        my $substr_off = $block_offset;
	my $len = $stat - $substr_off;

	#get entire feature
	$feat->{DATA} = substr($buffer->{DATA}, $substr_off, $len);
        $block_offset = $substr_off+$len; #end of current feature
	next if($meta && substr($feat->{DATA}, 0, 1) eq $meta); #skip comment line

	#get position part of feature entry (ignore the rest)
	$self->__fill_feature_ref($feat);

	#set file/buffer position info
	($feat->{FILE_OFFSET}, $feat->{BLOCK_OFFSET}) = $self->block_tell;
	$feat->{LENGTH} = $len; #length of total feature entry
	$buffer->{POS} += $feat->{LENGTH};

	$self->{FEATURE} = $feat;

	return $feat;
    }

    return;
}

#parse header from file
sub __parse_header {
    my $self = shift;
    my $IN = $self->__handle;
    my $header = $self->__header;

    #backup current position then move to start of file
    my $c_file_pos = $self->__file->{POS};
    my $c_real_pos = sysseek($IN, 0, SEEK_CUR);
    sysseek($IN, 0, SEEK_SET);

    #read header
    my @h_blocks;
    my $h_buffer = '';
    my $h_offset = 0;

    #skip line in header
    my $skip = $self->{IDX_SKIP};
    while($skip){
	#grow buffer
	my $eof;
	my $pos = index($h_buffer, "\n", $h_offset)+1;
	while(!$pos){
	    my $h_block = {};
	    $eof = ($self->__read_block($h_block)) ? 0 : 1;
	    last if($eof);
	    $h_block->{BUFFER_OFFSET} = length($h_buffer);
	    $header->{D_LENGTH} += $h_block->{D_LENGTH};
	    $self->__inflate(\ ($h_block->{DATA}), \$h_buffer);
	    $pos = index($h_buffer, "\n", $h_offset)+1;
	    push(@h_blocks, $h_block);
	}
	next if(!$pos && !$eof);
	
	#check
	die "ERROR: Malformed file\n" if(!$pos);
	$skip--;
	$h_offset = $pos;
    }

    #read comment lines
    my $meta = $self->{IDX_META};
    my $lead = substr($h_buffer, $h_offset, 1);
    while($meta && (!length($lead) || $lead eq $meta)){
	#grow buffer
	my $eof;
	my $pos = index($h_buffer, "\n", $h_offset)+1;
	while(!$pos){
	    my $h_block = {};
	    $eof = ($self->__read_block($h_block)) ? 0 : 1;
	    last if($eof);
	    $h_block->{BUFFER_OFFSET} = length($h_buffer);
	    $header->{D_LENGTH} += $h_block->{D_LENGTH};
	    $self->__inflate(\ ($h_block->{DATA}), \$h_buffer);
	    $pos = index($h_buffer, "\n", $h_offset)+1;
	    $lead ||= substr($h_buffer, $h_offset, 1);
	    push(@h_blocks, $h_block);
	    last if($lead ne $meta);
	}
	last if($lead ne $meta);
	next if(!$pos && !$eof);
	
	#check
	die "ERROR: Malformed file\n" if(!$pos);
	$h_offset = $pos;
	$lead = substr($h_buffer, $h_offset, 1);
    }

    #separate header from feature tail
    my $h_data = substr($h_buffer,0,$h_offset,'');
    $header->{TEXT} = $h_data;

    #identify real and virtual offset for end of header (start of features)
    if(!$self->{BGZIP}){ #uncompressed files
	$header->{TAIL_FILE_OFFSET} = length($h_data); #real offset
	$header->{TAIL_BLOCK_OFFSET} = 0; #virtual offset
	$header->{D_LENGTH} = length($h_data); #trim
	$h_buffer = ''; #no tail
    }
    else{
	my ($h_block) = grep {$h_offset < $_->{BUFFER_OFFSET} + $_->{I_LENGTH}} @h_blocks;
	my $b_off = $h_block->{BUFFER_OFFSET} || 0;
	$header->{TAIL_FILE_OFFSET} = $h_block->{FILE_OFFSET}; #real offset
	$header->{TAIL_BLOCK_OFFSET} = $h_offset-$b_off; #virtual offset
    }

    #deflate header (and tail)
    $self->__deflate(\ (substr($h_data, 0, 65280, '')), \ ($header->{DATA})) while(length($h_data));
    $self->__deflate(\ (substr($h_buffer, 0, 65280, '')), \ ($header->{TAIL})) while(length($h_buffer));

    #restore file positions
    $self->__file->{POS} = $c_file_pos;
    sysseek($IN, $c_real_pos, SEEK_SET);

    #parse reference info out of header
    $self->__fill_header_ref($header);

    #add unmapped
    $header->{NAME2SEQID}{'*'} = -1; #unmapped chr
    $header->{NAME2LENGTH}{'*'} = 0; #unmapped chr

    return $header;
}

sub __parse_idx {
    my $self = shift;
    my $idx_file = shift || $self->idx_file;

    #read in entire file
    my @tbi; #list representation
    my $raw; #raw data
    my $size = (stat($idx_file))[7];
    sysopen(my $IN, $idx_file, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $idx_file: $!";
    sysread($IN, $raw, $size);
    close($IN);

    #inflate tbi (it's BGZF compressed)
    my $data;
    while(length($raw)){
        #read BGZF header
	my $b_data = substr($raw, 0, 18, '');
        die "ERROR: Failure to read BGZF string\n" if(length($b_data) != 18);
	
	#validate header
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
    
    #check magic string
    my $str_off = 0; #offset in substring
    if(substr($data, $str_off, 4) ne "TBI\1"){
	die "ERROR: File does not appear to be a TBI index file\n"
    }
    $str_off += 4;
    
    #get number of reference sequences
    my $n_ref = unpack('l<', substr($data, $str_off, 4));
    $str_off += 4;

    #get formatting info
    my $format = unpack('l<', substr($data, $str_off, 4));
    $str_off += 4;
    my ($col_seq, $col_beg, $col_end) = unpack('l<l<l<', substr($data, $str_off, 12));
    $str_off += 12;
    my ($meta, $skip) = unpack('l<l<', substr($data, $str_off, 8));
    $str_off += 8;

    #get ref info
    my $l_nm = unpack('l<', substr($data, $str_off, 4));
    $str_off += 4;
    my @names = grep {length($_)} split(/\0/, substr($data, $str_off, $l_nm));
    my %name2tbid = map {$names[$_] => $_} (0..$#names);
    $str_off += $l_nm;

    #list of indices (n=n_ref)
    for(my $i = 0; $i < $n_ref; $i++){
	$tbi[$i] = {BINS => [], LINEAR => [], PSEUDO => undef}; #initialize index of each contig

	#get number of distinct bins (for the binning index)
	my $n_bin = unpack('l<', substr($data, $str_off, 4));
	$str_off += 4;
	
	#list of distinct bins (n=n_bin)
	my $tbi_b = $tbi[$i]{BINS};
	for(my $j = 0; $j < $n_bin; $j++){
	    #get distinct bin & number of chunks
	    my ($bin, $n_chunk) = unpack('L<l<', substr($data, $str_off, 8));
	    $str_off += 8;

	    #fill in pseudo-bin
	    if($bin == 37450){
		my $tbi_p = $tbi[$i]{PSEUDO} = []; #initialize chunk list for pseudo bin
		die "ERROR: Corrupt TBI index for pseudo bin\n" if($n_chunk != 2); #always 2

		#get (virtual) file offset of the start and end of placed unmapped reads
		my ($unmapped_beg, $unmapped_end) = unpack('Q<Q<', substr($data, $str_off, 16));
		$str_off += 16;

		#offsets are encoded as --> coffset<<16|uoffset
		my $coffset_beg = ($unmapped_beg >> 16);
		my $uoffset_beg = $unmapped_beg & 65535; #mask is (1<<16)-1
		my $coffset_end = ($unmapped_end >> 16);
		my $uoffset_end = $unmapped_end & 65535; #mask is (1<<16)-1

		push(@{$tbi_p}, [$coffset_beg, $uoffset_beg, $coffset_end, $uoffset_end]);

		#get number of mapped/unmapped read-segments for this reference
		my ($n_mapped, $n_unmapped) = unpack('Q<Q<', substr($data, $str_off, 16));
		$str_off += 16;

		push(@{$tbi_p}, [$n_mapped, $n_unmapped]);

		next;
	    }

	    #list of chunks (n=n_chunk)
	    $tbi_b->[$bin] = []; #initialize chunk list for each bin
	    for(my $k = 0; $k < $n_chunk; $k++){
		#get (virtual) file offset of the start and end of the chunk
		my ($chunk_beg, $chunk_end) = unpack('Q<Q<', substr($data, $str_off, 16));
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
	my $n_intv = unpack('l<', substr($data, $str_off, 4));
	$str_off += 4;
	
	#list of intervals (n=n_intv)
	my $tbi_l = $tbi[$i]{LINEAR};
	for(my $j = 0; $j < $n_intv; $j++){
	    #get (virtual) file offset of the first feature in the interval
	    my $ioffset = unpack('Q<', substr($data, $str_off, 8));
	    $str_off += 8;
	    
	    #interval offsets are encoded as --> coffset<<16|uoffset
	    my $i_coffset = ($ioffset >> 16);
	    my $i_uoffset = $ioffset & 65535; #mask is (1<<16)-1

	    $tbi_l->[$j] = [$i_coffset, $i_uoffset];
	}
    }
    
    #get number of unplaced unmapped reads (RNAME *)
    my $n_no_coor = unpack('Q<', substr($data, $str_off, 8));
    $str_off += 8;

    #hash representation
    my $tbi = $self->__idx;
    $tbi->{PATH} = $idx_file;
    $tbi->{ABS_PATH} = Cwd::abs_path($idx_file);
    $tbi->{LENGTH} = $size;
    $tbi->{TYPE} = $format;
    $tbi->{COL} = [$col_seq, $col_beg, $col_end];
    $tbi->{META} = $meta;
    $tbi->{SKIP} = $skip;
    $tbi->{INDICES} = \@tbi;
    $tbi->{N_NO_COOR} = $n_no_coor;
    $tbi->{NAMES} = \@names;
    $tbi->{NAME2TBID} = \%name2tbid;

    #add usful relation info from file header
    my $ref = $self->header->{REF};
    $tbi->{REF} = $ref;
    $tbi->{SEQID2TBID} = [map {$name2tbid{$self->seqid2name($_)}} (0..$#$ref)];
    $tbi->{TBID2SEQID} = [map {$self->name2seqid($names[$_])} (0..$#names)];

    return $n_ref;
}

sub __idx2string {
    my $self = shift;
    my $tbi = $self->tbi;

    #add magic
    my $data = "TBI\1";

    #number of reference sequences
    my $n_ref = scalar(@{$tbi->{NAMES}});
    $data .= pack('l<', $n_ref);

    #formatting info
    $data .= pack('l<', $tbi->{TYPE});
    $data .= pack('l<l<l<', @{$tbi->{COL}});
    $data .= pack('l<', $tbi->{META});
    $data .= pack('l<', $tbi->{SKIP});

    #ref info
    my $names = join("\0", @{$tbi->{NAMES}})."\0";
    my $l_nm = length($names);
    $data .= pack('l<', $l_nm);
    $data .= $names;

    #list of indices (n=n_ref)
    for(my $i = 0; $i < $n_ref; $i++){
	#number of distinct bins (for the binning index)
	my $tbi_b = $tbi->{INDICES}[$i]{BINS}; #real bins
	my $tbi_p = $tbi->{INDICES}[$i]{PSEUDO}; #pseudo bin
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
	my $tbi_l = $tbi->{INDICES}[$i]{LINEAR};
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
    my $n_no_coor = $tbi->{N_NO_COOR};
    $data .= pack('Q<', $n_no_coor) if(defined($n_no_coor));

    return $data;
}

sub seqid2name {
    my $self = shift;
    my $id = shift;
    my $header = $self->header;

    if(! exists($header->{REF}[$id])){
	if($self->{IDX}){
	    my $idx = $self->idx;
	    if(!exists($idx->{NAMES}[$id])){
		warn "WARNING: No contig exists with seqID $id\n";
		return undef;
	    }
	    
	    return $idx->{NAMES}[$id];
	}
	else{
	    warn "WARNING: No contig exists with seqID $id\n";
	    return undef;
	}
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

    if(!exists($header->{NAME2SEQID}{$name})){
	if($self->{IDX} && !$self->{IDX}{__HEADERFILL}){
	    my $idx = $self->idx;
	    my $names = $idx->{NAMES};
	    for(my $i = 0; $i < @$names; $i++){
		$header->{NAME2SEQID}{$names->[$i]} ||= $i;
	    }
	    $self->{IDX}{__HEADERFILL} = 1; #flag to only use IDX to fill header once

	    if(!exists($header->{NAME2SEQID}{$name})){
		warn "WARNING: No contig exists with name $name\n";
		return undef;
	    }
	    
	    return $header->{NAME2SEQID}{$name};
	}
	else{
	    warn "WARNING: No contig exists with name $name\n";
	    return undef;
	}
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

    #not compreseed
    if(!$self->{BGZIP}){
	$$out_buffer .= $$in_buffer;
	$$in_buffer = '';
	return;
    }

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

    #not compreseed
    if(!$self->{BGZIP}){
	$$out_buffer .= $$in_buffer;
	$$in_buffer = '';
        return;
    }

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

    #not compreseed
    if(!$self->{BGZIP}){
	my $stat = sysread($IN, $block->{DATA}, 65536, 0);
	return undef if($stat == 0); #end of file
	die "ERROR: Failure to read file stream\n" if($stat < 0);

	#fill in block values
	my $file = $self->__file;
	$block->{FILE_OFFSET} = $file->{POS};
	$block->{D_LENGTH} = $stat;
	$block->{I_LENGTH} = $stat;
	$file->{POS} += $block->{D_LENGTH}; #move file position

	return 1;
    }

    #read BGZF header
    my $stat = sysread($IN, $block->{DATA}, 18, 0);
    return undef if($stat == 0); #end of file
    die "ERROR: Failure to read file stream\n" if($stat != 18);

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
    my $file = $self->__file;
    my $buffer = $self->__buffer;

    #go to position
    my $pos = $self->__file->{POS};
    my $size = $self->__file->{LENGTH}; #file size

    #start of file
    if($pos == 0){
	$buffer->{UPSTREAM} = "\n";
	$file->{POS} = sysseek($IN, $pos, 0);
	return $file->{POS};
    }
    elsif(!$self->{BGZIP}){ #not compressed
	sysseek($IN, $pos-1, SEEK_SET); #one character upstream
	my $stat = sysread($IN, $buffer->{UPSTREAM}, 1, 0);
	die "ERROR: Failure to read file: $!\n" if($stat < 0);
	$file->{POS} = sysseek($IN, $pos, 0);
	return $file->{POS};
    }

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
        if($offset == -1){
	    $file->{POS} = sysseek($IN, $size, 0);
	    return $file->{POS};
	}

        #match second static key
        my $offset2 = index($data, $key2, $offset);
        if($offset2 == -1){
            $file->{POS} = sysseek($IN, $size, 0);
            return $file->{POS};
	}

        #second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
        if($offset2-$offset == 10){
            last;
        }
        else{
            $offset += 4; #jump over key1 match and try again
        }
    }
    $pos += $offset; #adjust file position with offset

    #now find block imediately prior to this one
    $offset = 0;
    my $alt = ($pos-65536 > 0) ? $pos-65536 : 0; #max block backwards
    sysseek($IN, $alt, SEEK_SET);
    sysread($IN, $data, $pos-$alt, 0); #everything up to current block
    while(1){
        #match first static key
        $offset = index($data, $key1, $offset);
        if($offset == -1){
            $file->{POS} = sysseek($IN, $size, 0);
            return $file->{POS};
	}

        #match second static key
        my $offset2 = index($data, $key2, $offset);
        if($offset2 == -1){
            $file->{POS} = sysseek($IN, $size, 0);
            return $file->{POS};
	}

        #second match must be 10 bytes from first offset (4 byte key1 and 6 byte gap)
        if($offset2-$offset == 10){
	    last;
        }
        else{
            $offset += 4; #jump over key1 match and try again
        }
    }
    substr($data, 0, $offset, ''); #trim

    #inflate one block at a time
    my $upstream = '';
    while(length($data)){
        #read BGZF header
        my $b_data = substr($data, 0, 18, '');
        die "ERROR: Failure to read BGZF string\n" if(length($b_data) != 18);

        #validate header
        my ($id1,$id2,$cm,$flg,$mtime,$xfl,$os,
            $xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $b_data);
        die "ERROR: Does not appear to be a BGZF file\n"
            if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
               $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);

        #read remainder compression block and footer
        my $c_data_size = $bsize-$xlen-19; #compression block size
        $b_data .= substr($data, 0, $c_data_size + 8, ''); #the +8 is for the footer
        die "ERROR: Could not read compression block\n"
            if(length($b_data) != 18 + $c_data_size + 8);

        #inflate and concatenate to output
        $self->__inflate(\$b_data, \$upstream);
    }

    #get last upstream character
    $buffer->{UPSTREAM} = substr($upstream, -1, 1);
    $buffer->{UPSTREAM} = "\n" if(!length($buffer->{UPSTREAM}));

    #return
    $file->{POS} = sysseek($IN, $pos, 0);
    return $file->{POS};
}

#walks across uncompressed block looking for feature header
sub __adjust_to_next_feature {
    my $self = shift;
    my $buffer = $self->__buffer;
    my $header = $self->header;

    #reset value
    $self->{FEATURE} = undef if($self->{FEATURE});

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

    #adjust for header and missing data
    my ($f_off, $b_off) = $self->block_tell;
    my ($tf_off, $tb_off) = ($header->{TAIL_FILE_OFFSET}, $header->{TAIL_BLOCK_OFFSET});
    if($f_off < $tf_off || ($f_off == $tf_off && $b_off <= $tb_off)){ #jump to end of header
	$self->block_seek($tf_off, $tb_off);
	return 0;
    }
    elsif($buffer->{POS} == 0 && $buffer->{UPSTREAM} eq "\n"){ #previous block ends on endline
        return 0; #seems to be a valid line
    }

    #search for downstream endline
    while(1){
	my $pos = $buffer->{POS};
	if(my $next = index($buffer->{DATA}, "\n", $pos)+1){
	    $buffer->{POS} = $next;
	    return $next;
	}
	else{
	    return unless($self->__read_block);
	    $self->__inflate_to_buffer;
	}
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
	    my $upstream = substr($buffer->{DATA}, 0, $block->{I_LENGTH}, ''); #remove text from buffer
	    $buffer->{POS} -= $block->{I_LENGTH};
	    $buffer->{LENGTH} -= $block->{I_LENGTH};
	    $buffer->{UPSTREAM} = substr($upstream, -1, 1) if(length($upstream));
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

#------------------------------------------------------------------------
#--------------                 FUNCTIONS                  --------------
#------------------------------------------------------------------------

#end of file signature block
sub eof_block {
    return pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000');
}

sub pos2bin {
    return (( ((1<<15)-1)/7 + ($_[0]>>14) ) % 65536);
}

sub pos2bins {
    return (0,
	    (   1 + ($_[0]>>26)) % 65536,
	    (   9 + ($_[0]>>23)) % 65536,
	    (  73 + ($_[0]>>20)) % 65536,
	    ( 585 + ($_[0]>>17)) % 65536,
	    (4681 + ($_[0]>>14)) % 65536);
}

sub reg2bin {
    my $beg = shift;
    my $end = shift;

    $end-- if($beg != $end);
    return (( ((1<<15)-1)/7 + ($beg>>14) ) % 65536) if($beg>>14 == $end>>14);
    return (( ((1<<12)-1)/7 + ($beg>>17) ) % 65536) if($beg>>17 == $end>>17);
    return (( ((1<<9)-1)/7  + ($beg>>20) ) % 65536) if($beg>>20 == $end>>20);
    return (( ((1<<6)-1)/7  + ($beg>>23) ) % 65536) if($beg>>23 == $end>>23);
    return (( ((1<<3)-1)/7  + ($beg>>26) ) % 65536) if($beg>>26 == $end>>26);
    return 0;
}

sub reg2bins {
    my $beg = shift;
    my $end = shift;

    $end-- if($beg != $end);
    my $i = 1;
    my @list = 0;
    for(my $k =    1 + ($beg>>26); $k <=    1 + ($end>>26); $k++){$list[$i++] = $k % 65536}
    for(my $k =    9 + ($beg>>23); $k <=    9 + ($end>>23); $k++){$list[$i++] = $k % 65536}
    for(my $k =   73 + ($beg>>20); $k <=   73 + ($end>>20); $k++){$list[$i++] = $k % 65536}
    for(my $k =  585 + ($beg>>17); $k <=  585 + ($end>>17); $k++){$list[$i++] = $k % 65536}
    for(my $k = 4681 + ($beg>>14); $k <= 4681 + ($end>>14); $k++){$list[$i++] = $k % 65536}

    return \@list;
}

sub _cap_hash {
    my $rhash = shift;
    my %hash = map {
        my $k = $_;
	my $v = $rhash->{$k};
        $k =~ tr/a-z/A-Z/;
	$k => $v;
    } keys(%{$rhash});
    return \%hash;
}

1;
