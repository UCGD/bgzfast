=head1 NAME

BGZFast::bgzf_utility - Convenient functions BGZF file manipulation

=head1 SYNOPSIS

  use BGZFast::bgzf_utility;

  # Find first block in a string
  my $offset = BGZFast::bgzf_utility::find_block_offset($str, 0);

  # Find last block in a string
  $offset = BGZFast::bgzf_utility::rfind_block_offset($str, length($string));

  #find full blocks
  my ($off, $len) = BGZFast::bgzf_utility::find_full_block($str, 0); #first
  ($off, $len) = BGZFast::bgzf_utility::rfind_full_block($str, length($string)); #last

  # Get info for a block
  my ($blen, $ilen, $crc32) = BGZFast::bgzf_utility::block_stat($offset);

  # Compress/Decompress a string
  my $unzipped = BGZFast::bgzf_utility::unzip($str);
  my $zipped = BGZFast::bgzf_utility::zip($unzipped);

=head1 DESCRIPTION

BGZFast::bgzf_utility provides convenient functions for basic manipulation of
BGZF file structure and content.

=head1 BASICS OF BGZF HEADER/FOOTER UNPACKING AND VALIDATION

BGZFast::bgzf_utility manipulates and validates BGZF blocks using specific
assumptions based on the format specification. The basic logic of these
assumptions is simple but may be unclear when reading the source code
directly. The following non-optimized pseudo-code demostrates the same basic
logic in a less obfuscated and hopfully easy to understand way.

  # In BGZF format data is compressed into a series of indepenent blocks.
  # The header for each block is 18 bytes long and can be unpacked like so
  my $header = substr($block, 0, 18); #first 18 bytes of block
  my ($id1,$id2,$cm,$flg, $mtime,$xfl,$os,
      $xlen, $si1,$si2,$slen,$bsize) = unpack('CCCCVCCvCCvv', $header);

  # Most values in block header are static and can be validated like so
  die "ERROR: Does not appear to be a BGZF file\n"
     if($id1 != 31 || $id2 != 139 || $cm != 8 || $flg != 4 ||
        $xlen != 6 || $si1 != 66 || $si2 != 67 || $slen != 2);

  # The value of $bsize cannot be larger than 65536 bytes or 64kb
  die "ERROR: Bsize is outside of expected range\n" if($bsize > 65536);

  # Note $bsize reported in the header is always 1 shorter than true value
  my $true_size = $bsize + 1;
  die "ERROR: Calculated block size does not match observed size\n"
     if($true_size != length($block));
 
  # Calculate the size of the compressed data section following the header
  my $cdata_size = $bsize - $xlen - 19; #Note $xlen is always 6
  my $cdata = substr($block, 18, $cdata_size); #offset of 18 to skip header
 
  # Block footer is 8 bytes long and can be unopacked like so
  my $footer = substr($block, 18 + $cdata_size, 8); #follows compressed section
  my ($crc32, $isize) = unpack('L<L<', $footer);

  # The isize or inflated data size should be smaller than bsize to
  # ensure that if it grows with compression (can happen), then it does
  # not force bsize to exceed its maximum. Samtools does this by setting
  # isize to a maximum of 65536 minus 256 or 65280.
  die "ERROR: Isize is outside of expected range\n" if($isize > 65280);

  # The isize value and crc32 value can be validated after decompression
  my $data = Compress::Zlib::uncompress($cdata);
  die "ERROR: Inflated data does not match expected size\n"
     if(length($data) != $isize);
  die "ERROR: Inflated data CRC32 not match expected value\n"
     if(Compress::Zlib::crc32($data) != $crc32);

=head1 SEE ALSO

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
#----                     BGZFast::bgzf_utility                      ----
#------------------------------------------------------------------------
package BGZFast::bgzf_utility;

use strict;
use warnings;
use Compress::Zlib;
use Fcntl qw(:DEFAULT :seek); #SEEK_SET, SEEK_CUR, SEEK_END
use constant SEEK_HOLE => 4; #SEEK_HOLE

our $VERSION = '0.01';

#constants
use constant KEY1 => pack('H8', '1f8b0804'); #4 byte static header key
use constant KEY2 => pack('H12', '060042430200'); #key1, 6 non-static bytes, then key2 (6 bytes)
use constant EOF_BLOCK  => pack('H56', '1f8b08040000000000ff0600424302001b0003000000000000000000'); #EOF block
use constant E_HEADER  => pack('H36', '1f8b08040000000000ff0600424302000000'); #empty BGZF header
use constant MAX_BSIZE => 65536; #max BSIZE is 64kb or 2^16 bytes (18 header + CDATA_SIZE + 8 footer + 1 extra)
use constant MAX_ISIZE => 65536-256; #ensures that BSIZE is never becomes too big for "bigger when compressed" blocks

#------------------------------------------------------------------------
=head2 validate_block_header

 Title   : validate_block_header
 Usage   : BGZFast::bgzf_utility::validate_block_header($string);
 Function: Validate the BGZF header in the given string
 Returns : True on success, 0 on failure, and undef if string too short
 Args    : String or string reference
           Optional offset in string

=cut
sub validate_block_header {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $o = $_[1] || 0;
    return undef if(length($$ref)-$o < 16);
    return (substr($$ref, $o, 4) eq KEY1 && substr($$ref, $o+10, 6) eq KEY2) ? 1 : 0;
}
#------------------------------------------------------------------------
=head2 block_stat

 Title   : block_stat
 Usage   : BGZFast::bgzf_utility::block_stat($string);
 Function: Get information about a BGZF block from its header and footer 
 Returns : Block length, inflated data length, and CRC32 of inflated data
           (will return undef on partial block).
 Args    : String or string reference
           Optional offset in string

=cut
sub block_stat {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $offset = $_[1] || 0;
    return undef if(length($$ref) < 28); #not enough for full block
    die "ERROR: Invalid block header\n" if(!validate_block_header($ref, $offset));
    my $bsize = scalar(unpack('v', substr($$ref, $offset+16, 2))) + 1; #+1 because bsize is always too small
    die "ERROR: Invalid BSIZE. Block may be corrupt.\n" if($bsize > MAX_BSIZE);
    return undef if(length($$ref) < $offset+$bsize); #not enough for full block
    my ($crc32, $isize) = unpack('L<L<', substr($$ref, $offset+$bsize-8, 8)); #footer
    warn "WARNING: ISIZE is too large. Block may be corrupt.\n" if($isize > MAX_ISIZE); #no die since technically legal
    return ($bsize, $isize, $crc32);
}
#------------------------------------------------------------------------
=head2 find_block_offset

 Title   : find_block_offset
 Usage   : BGZFast::bgzf_utility::find_block_offset($string);
 Function: Scan for start of BGZF block in string
 Returns : Offset immediately before block start (-1 if not found)
 Args    : String or string reference
           Optional offset in string

=cut
sub find_block_offset {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $offset = $_[1] || 0;

    $offset = 0 if($offset < 0);

    #find block keys in data
    my $max = length($$ref)-16; #need at least 16 bytes for a match
    while(1){
	return -1 if($offset > $max);

        #match first static key
        $offset = index($$ref, KEY1, $offset);
        return -1 if($offset == -1);

        #match second static key
	if(substr($$ref, $offset+10, 6) eq KEY2){ #key1 (4 bytes), 6 other bytes, key2 (6 bytes)
	    last;
        }
        else{
            $offset += 4; #jump over key1 match window and try again
        }
    }

    return $offset;
}
#------------------------------------------------------------------------
=head2 find_full_block

 Title   : find_full_block
 Usage   : BGZFast::bgzf_utility::find_full_block($string);
 Function: Scan for the first full BGZF block in a string
 Returns : Offset of block start and block length (returns -1 if not found)
 Args    : String or string reference
           Optional offset in string

=cut
sub find_full_block {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    return -1 if(length($$ref) < 28); #minimum block is 28 bytes

    my $offset = find_block_offset(@_);
    return -1 if($offset == -1);
    return -1 if($offset > length($ref)-28); #minimum block is 28 bytes

    #validate full block contained
    my $bsize = scalar(unpack('v', substr($$ref, $offset+16, 2))) + 1; #+1 because bsize is always too small
    return -1 if($offset+$bsize > length($ref)); 
    return wantarray ? ($offset, $bsize) : $offset;
}
#------------------------------------------------------------------------
=head2 rfind_block_offset

 Title   : rfind_block_offset
 Usage   : BGZFast::bgzf_utility::rfind_block_offset($string);
 Function: Scan backwards for start of BGZF block in string
 Returns : Offset immediately before block start (-1 if not found)
 Args    : String or string reference
           Optional offset in string

=cut
sub rfind_block_offset {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $offset = $_[1] || 0;

    $offset = length($$ref)-16 if($offset > length($$ref)-16); #need at least 16 bytes to match

    #find block keys in data
    while(1){
	return -1 if($offset < 0); #min

        #match first static key
        $offset = rindex($$ref, KEY1, $offset);
        return -1 if($offset == -1);

        #match second static key
	if(substr($$ref, $offset+10, 6) eq KEY2){ #key1 (4 bytes), 6 other bytes, key2 (6 bytes)
            last;
        }
        else{
            $offset -= 4; #jump over key1 match window and try again
        }
    }

    return $offset;
}
#------------------------------------------------------------------------
=head2 rfind_full_block

 Title   : rfind_full_block
 Usage   : BGZFast::bgzf_utility::rfind_full_block($string);
 Function: Scan backwards for the first full BGZF block in a string
 Returns : Offset of block start and block length (returns -1 if not found)
 Args    : String or string reference
           Optional offset in string

=cut
sub rfind_full_block {
    my $ref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    return -1 if(length($$ref) < 28); #minimum block is 28 bytes

    my $offset = $_[1] || 0;
    $offset = length($$ref)-28 if($offset > length($$ref)-28);
    $offset = rfind_block_offset($ref, $offset);
    return -1 if($offset == -1);

    #validate full block contained
    my $bsize = scalar(unpack('v', substr($$ref, $offset+16, 2))) + 1; #+1 because bsize is always too big
    $offset = rfind_block_offset($ref, $offset-4) if($offset+$bsize > length($ref));
    return wantarray ? ($offset, $bsize) : $offset;
}
#------------------------------------------------------------------------
=head2 zip

 Title   : zip
 Usage   : $zipped = BGZFast::bgzf_utility::zip($string);
           BGZFast::bgzf_utility::zip($string, $zipped, 6);
 Function: Deflate a string to BGZF format
 Returns : String of deflated data
 Args    : String or string reference of data to deflate
           Optional string/ref to deflate to (returns new otherwise)
           Optional desired compression level (defaults to 6)

 Note that this method is destuctive to input string but concatenates to
 the output string.

=cut
sub zip {
    my $iref  = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $oref  = $_[1];
    my $level = $_[2];
    
    if(!$oref){
        my $string = '';
        $oref = \$string;
    }
    elsif(!ref($oref)){
	$oref = \ ($_[1]);
    }

    return $$oref if(!length($$iref));

    #compress data from string to string
    while(length($$iref)){
        __deflate(\ (substr($$iref, 0, MAX_ISIZE)), $oref, $level);
    }

    return $$oref;
}
#------------------------------------------------------------------------
=head2 unzip

 Title   : unzip, unzip_partial, unzip_full
 Usage   : $string = BGZFast::bgzf_utility::unzip($zipped);
           BGZFast::bgzf_utility::unzip($zipped, $string, 1); #partial fail
           BGZFast::bgzf_utility::unzip_partial($zipped, $string);
           BGZFast::bgzf_utility::unzip_full($zipped, $string); #partial fail
 Function: Deflate a string to BGZF format
 Returns : String of deflated data
 Args    : String or string reference of data to deflate
           Optional string/ref to deflate to (returns new otherwise)
           Optional flag to induce failure on partial blocks (default false)
        
 Note that this method is destuctive to input string but concatenates to
 the output string. Will leave trailing partial blocks in input string if
 failure flag not set.

=cut
sub unzip_full {return unzip($_[0], $_[1], 1)} #fails on partial
sub unzip_partial {return unzip($_[0], $_[1], 0)} #does not fail on partial
sub unzip {
    my $iref = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $oref = $_[1];
    my $fail = $_[2] || 0;  #whether to fail on partial blocks

    if(!$oref){
        my $string = '';
        $oref = \$string;
    }
    elsif(!ref($oref)){
	$oref = \ ($_[1]);
    }

    return $$oref if(!length($$iref));
    
    #decompress data from string to string
    my $bsize; #predeclare
    while(length($$iref) >= 28){
        #validate BGZF header
	die "ERROR: Does not appear to be a BGZF file\n"
	    if(substr($$iref, 0, 4) ne KEY1 || substr($$iref, 10, 6) ne KEY2);

	#compression block size (bsize-xlen-19)
	$bsize = scalar(unpack('v', substr($$iref, 16, 2))) + 1; #+1 because BSIZE is always too small 
	die "ERROR: Invalid BSIZE. Block may be corrupt.\n" if($bsize > MAX_BSIZE);
        last if(length($$iref) < $bsize); #check if string contains full block

        #inflate and concatenate to output
        __inflate(\ (substr($$iref, 0, $bsize)), $oref);
    }

    #fail on trailing block if indicated
    die "ERROR: Could not read compression block\n"
	if($fail && length($$iref));

    return $$oref;
}
#------------------------------------------------------------------------
#inflate a single compression block (very limited error checking)
sub __inflate {
    my $in_buffer  = $_[0]; #reference
    my $out_buffer = $_[1]; #reference

    return if(!length($$in_buffer));

    my ($i_obj, $stat) = inflateInit(-WindowBits => -15, -Bufsize => 131072);
    die "ERROR: Failed to create zlib inflation object with status: $stat\n" if($stat);

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
#------------------------------------------------------------------------
#deflate a single compresssion block
sub __deflate {
    my $in_buffer = $_[0]; #reference
    my $out_buffer = $_[1]; #reference
    my $level = $_[2];

    $level = 6 if(!defined($level) || $level < 0);
    
    return EOF_BLOCK if(!length($$in_buffer));

    my ($d_obj, $stat) = deflateInit(-WindowBits => -15, -Bufsize => 131072, -Level => $level);
    die "ERROR: Failed to create zlib deflation object with status: $stat\n" if($stat);

    my $offset = length($$out_buffer);
    $$out_buffer .= E_HEADER; #header with BSIZE=0
    $$out_buffer .= scalar($d_obj->deflate($in_buffer)).scalar($d_obj->flush); #compressed data
    $$out_buffer .= pack('V', crc32($in_buffer)); #CRC32
    $$out_buffer .= pack('V', $d_obj->total_in()); #ISIZE
    substr($$out_buffer, $offset+16, 2, pack('v', $d_obj->total_out+6+19)); #set final BSIZE
    $$in_buffer = ''; #destroy input buffer

    return;
}
#------------------------------------------------------------------------
=head2 sparse_check

 Title   : sparse_check
 Usage   : BGZFast::bgzf_utility::sparse_check($file) or die "corrupt file";
 Function: Checks for sparse file holes which will be null padded. 
 Returns : 1 if a hole is found and 0 otherwise (undef if not supported)
 Args    : File path

=cut
sub sparse_check {
    my $file = $_[0];

    my $size = (stat($file))[7];
    return 0 if(!$size);

    #null padding validation
    sysopen(my $IN, $file, O_RDONLY|O_BINARY) or die "ERROR: Could not open file $file: $!";
    my $hole = sysseek($IN, 0, SEEK_HOLE); #4 is SEEK_HOLE
    close($IN);

    #OS does not support SEEK_HOLE
    if(!defined($hole)){
	#warn "WARNING: OS does not support the SEEK_HOLE flag\n";
	return undef;
    }

    return (defined($hole) && $size != $hole) ? 1 : 0;
}
#------------------------------------------------------------------------
=head2 eof_block

 Title   : eof_block
 Usage   : $eof = BGZFast::bgzf_utility::eof_block();
 Function: Get a standard BGZF EOF block
 Returns : String of EOF block
 Args    : None

=cut
sub eof_block {
    return EOF_BLOCK;
}
#------------------------------------------------------------------------
=head2 max_bsize

 Title   : max_bsize
 Usage   : $eof = BGZFast::bgzf_utility::max_bsize();
 Function: Get the maximum block size as defined for BGZF format
 Returns : 65536 (== 2**16 == 64kb)
 Args    : None

=cut
sub max_bsize {
    return MAX_BSIZE;
}
#------------------------------------------------------------------------
=head2 max_isize

 Title   : max_isize
 Usage   : $eof = BGZFast::bgzf_utility::max_isize();
 Function: Get the maximum size of data to compress for a single block
           Must be smaller than MAX_BSIZE because it is possible for
           blocks to grow with compression (samtools uses MAX_BSIZE-256)
 Returns : 65280 (MAX_BSIZE - 256)
 Args    : None

=cut
sub max_isize {
    return MAX_ISIZE;
}
#------------------------------------------------------------------------
1;
