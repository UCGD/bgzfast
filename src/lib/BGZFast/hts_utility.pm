=head1 NAME

BGZFast::hts_utility - Convenient functions for using HTS bins

=head1 SYNOPSIS

  use BGZFast::hts_utility;

  # Get smallest bin completely encapsulating region
  my $bin = BGZFast::hts_utility::reg2bin($start, $end);

  # Get all regions overlapping a region
  my @bins = BGZFast::hts_utility::reg2bins($start, $end);

  # Get region for a bin
  ($start, $end) = BGZFast::hts_utility::bin2reg($bin);

  # Get immediate parent of a bin
  my $parent = BGZFast::hts_utility::bin2paren($bin);

=head1 DESCRIPTION

BGZFast::hts_utility provides convenient functions for using HTS bins.

=head1 BASICS OF HTS BIN USAGE AND CALCULATION

BGZFast::hts_utility calculates bins and relationships among bins.
The binning scheme used is described in the BAM Index format documentation.
Positions and regions used for binning are zero based with a closed-open
interval [start,end). This means the start is included but the end is not.
Basic theory is as follows:

  #Bin levels are based on bit shifts (max 29 and min 14)
  29, 26, 23, 20, 17, 14

  #Bin count shifts for each level are calculated as so
  shift = ((1<< 0)-1)/7 =    0
  shift = ((1<< 3)-1)/7 =    1
  shift = ((1<< 6)-1)/7 =    9
  shift = ((1<< 9)-1)/7 =   73
  shift = ((1<<12)-1)/7 =  585
  shift = ((1<<15)-1)/7 = 4681

  #Binning at each level is done by bit shifting then count shifting
  bin = (pos>>bit) + shift
  bin = (pos>>29)  +    0
  bin = (pos>>26)  +    1
  bin = (pos>>23)  +    9
  bin = (pos>>20)  +   73
  bin = (pos>>17)  +  585
  bin = (pos>>14)  + 4681

  #Note that linear index bins are 2^14 bytes wide and never shifted
  lin = (pos>>14)

=head1 SEE ALSO

L<BGZFast::tbi_manipulator>

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
#----                      BGZFast::hts_utility                      ----
#------------------------------------------------------------------------
package BGZFast::hts_utility;

use strict;
use warnings;

use constant MIN_BINSIZE => 16384;     #2**14
use constant MAX_BINSIZE => 536870912; #2**29
use constant PSEUDO_BIN  => 37450;     #Max possible bin + 2

our $VERSION = '0.01';


#------------------------------------------------------------------------
#----                            GENERAL                             ----
#------------------------------------------------------------------------
=head2 pseudo_bin

 Title   : pseudo_bin
 Usage   : my $pseudo = BGZFast::hts_utility::pseudo_bin();
 Function: Get the pseudo bin number. This occurs two bins beyond the
           maximum possible bin number.
 Returns : Pseudo bin number
 Args    : None

=cut
sub pseudo_bin {
    return PSEUDO_BIN;
}
#------------------------------------------------------------------------
=head2 min_binsize

 Title   : min_binsize
 Usage   : my $min = BGZFast::hts_utility::min_binsize();
 Function: Get the minimum size of a standard bin. Corresponds to size
           of bin 4681 and above.
 Returns : Minimum size of a standard bin
 Args    : None

=cut
sub min_binsize {
    return MIN_BINSIZE;
}
#------------------------------------------------------------------------
=head2 max_binsize

 Title   : max_binsize
 Usage   : my $max = BGZFast::hts_utility::max_binsize();
 Function: Get the maximum size of a standard bin. Corresponds to size
           of bin 0.
 Returns : Maximum size of a standard bin
 Args    : None

=cut
sub max_binsize {
    return MAX_BINSIZE;
}
#------------------------------------------------------------------------
=head2 bin2bit

 Title   : bin2bit
 Usage   : my $bit = BGZFast::hts_utility::bin2bit($bin);
 Function: Get the bit shift used to calculate bin size/width (log2 of
           bin size/width)
 Returns : Bit shift number
 Args    : Bin number

=cut
sub bin2bit {
    return 14 if($_[0] >= 4681); #2**14
    return 17 if($_[0] >=  585); #2**17
    return 20 if($_[0] >=   73); #2**20
    return 23 if($_[0] >=    9); #2**23
    return 26 if($_[0] >=    1); #2**26
    return 29 if($_[0] ==    0); #2**29
}
#------------------------------------------------------------------------
=head2 bin2size

 Title   : bin2size
 Usage   : my $bsize = BGZFast::hts_utility::bin2size($bin);
 Function: Get the size/width of a goiven bin. For example bin 9 has a
           width of 8Mb and bin 4681 has a width of 16kb
 Returns : The bin size/width
 Args    : Bin number

=cut
sub bin2size {
    return     16384 if($_[0] >= 4681); #2**14
    return    131072 if($_[0] >=  585); #2**17
    return   1048576 if($_[0] >=   73); #2**20
    return   8388608 if($_[0] >=    9); #2**23
    return  67108864 if($_[0] >=    1); #2**26
    return 536870912 if($_[0] ==    0); #2**29
}
#------------------------------------------------------------------------
=head2 bin2shift

 Title   : bin2shift
 Usage   : my $bshift = BGZFast::hts_utility::bin2shift($bin);
 Function: Get the bin shift added for a given bin. For example 8Mb bins
           start at 9 and all 16kb bins start at 4681  
 Returns : The bin count shift added to the given bin
 Args    : Bin number

=cut
sub bin2shift {
    return 4681 if($_[0] >= 4681); #((1<<15)-1)/7
    return  585 if($_[0] >=  585); #((1<<12)-1)/7
    return   73 if($_[0] >=   73); #((1<< 9)-1)/7
    return    9 if($_[0] >=    9); #((1<< 6)-1)/7
    return    1 if($_[0] >=    1); #((1<< 3)-1)/7
    return    0 if($_[0] ==    0); #((1<< 0)-1)/7
}

#------------------------------------------------------------------------
#----                          LINEAR INDEX                          ----
#------------------------------------------------------------------------
=head2 pos2lin

 Title   : pos2lin
 Usage   : my $lin = BGZFast::hts_utility::pos2lin($pos);
 Function: Get linear index bin containing position.
 Returns : Linear index bin number
 Args    : Position

=cut
sub pos2lin {
    return ($_[0]>>14) % 65536;
}
#------------------------------------------------------------------------
=head2 reg2lins

 Title   : reg2lins
 Usage   : my @lins = BGZFast::hts_utility::reg2lins($start, $end);
 Function: Get all linear index bins overlapping a region
 Returns : Linear index bin numbers
 Args    : Start and end

=cut
sub reg2lins {
    $_[1]--;
    my $i = 0;
    my @list;
    for(my $k = ($_[0]>>14); $k <= ($_[1]>>14); $k++){$list[$i++] = $k % 65536}

    return \@list;
}

#------------------------------------------------------------------------
#----                           BIN INDEX                            ----
#------------------------------------------------------------------------
=head2 pos2bin

 Title   : pos2bin
 Usage   : my $bin = BGZFast::hts_utility::pos2bin($pos);
 Function: Get smallest standard bin containing position. Will always be
           a value greater than 4680.
 Returns : Bin number
 Args    : Position

=cut
sub pos2bin {
    return ( 4681 + ($_[0]>>14) ) % 65536;
}
#------------------------------------------------------------------------
=head2 pos2bins

 Title   : pos2bins
 Usage   : my @bins = BGZFast::hts_utility::pos2bins($pos);
 Function: Get all standard bins overlapping a position.
 Returns : Bin numbers
 Args    : Position

=cut
sub pos2bins {
    return (0,
            (   1 + ($_[0]>>26)) % 65536,
            (   9 + ($_[0]>>23)) % 65536,
            (  73 + ($_[0]>>20)) % 65536,
            ( 585 + ($_[0]>>17)) % 65536,
            (4681 + ($_[0]>>14)) % 65536);
}
#------------------------------------------------------------------------
=head2 reg2bin

 Title   : reg2bin
 Usage   : my $bin = BGZFast::hts_utility::reg2bin($start, $end);
 Function: Get smallest standard bin fully containing a region
 Returns : Bin number
 Args    : Start and end

=cut
sub reg2bin {
    $_[1]-- if($_[0] != $_[1]);
    return (( 4681 + ($_[0]>>14) ) % 65536) if($_[0]>>14 == $_[1]>>14);
    return ((  585 + ($_[0]>>17) ) % 65536) if($_[0]>>17 == $_[1]>>17);
    return ((   73 + ($_[0]>>20) ) % 65536) if($_[0]>>20 == $_[1]>>20);
    return ((    9 + ($_[0]>>23) ) % 65536) if($_[0]>>23 == $_[1]>>23);
    return ((    1 + ($_[0]>>26) ) % 65536) if($_[0]>>26 == $_[1]>>26);
    return 0;
}
#------------------------------------------------------------------------
=head2 reg2bins

 Title   : reg2bins
 Usage   : my @bins = BGZFast::hts_utility::reg2bins($start, $end);
 Function: Get all standard bins overlapping a region
 Returns : Bin numbers
 Args    : Start and end

=cut
sub reg2bins {
    $_[1]-- if($_[0] != $_[1]);
    my $i = 1;
    my @list = 0;
    for(my $k =    1 + ($_[0]>>26); $k <=    1 + ($_[1]>>26); $k++){$list[$i++] = $k % 65536}
    for(my $k =    9 + ($_[0]>>23); $k <=    9 + ($_[1]>>23); $k++){$list[$i++] = $k % 65536}
    for(my $k =   73 + ($_[0]>>20); $k <=   73 + ($_[1]>>20); $k++){$list[$i++] = $k % 65536}
    for(my $k =  585 + ($_[0]>>17); $k <=  585 + ($_[1]>>17); $k++){$list[$i++] = $k % 65536}
    for(my $k = 4681 + ($_[0]>>14); $k <= 4681 + ($_[1]>>14); $k++){$list[$i++] = $k % 65536}

    return \@list;
}

#------------------------------------------------------------------------
#----                       USEFUL CONVERSIONS                       ----
#------------------------------------------------------------------------
=head2 lin2bin

 Title   : lin2bin
 Usage   : my $bin = BGZFast::hts_utility::lin2bin($lin);
 Function: Get the standard bin equivilent to linear index bin.
 Returns : Bin number
 Args    : Linear index bin

=cut
sub lin2bin {
    return $_[0] + 4681;
}
#------------------------------------------------------------------------
=head2 bin2lin

 Title   : bin2lin
 Usage   : my $lin = BGZFast::hts_utility::bin2lin($bin);
 Function: Get the linear index bin equivilent to standard bin. Only
           works for bin numbers higher than 4680.
 Returns : Linear index bin
 Args    : Bin number

=cut
sub bin2lin {
    return $_[0] - 4681 if($_[0] >= 4681);
    die "ERROR: No direct linear index conversion. Try bin2lins instead\n";;
}
#------------------------------------------------------------------------
=head2 bin2lins

 Title   : bin2lins
 Usage   : my @lins = BGZFast::hts_utility::bin2lins($bin);
 Function: Get all linear index bins overlapping standard bins
 Returns : Linear index bins
 Args    : Bin number

=cut
sub bin2lins {
    return reg2lins(bin2reg($_[0]));
}
#------------------------------------------------------------------------
=head2 bin2beg

 Title   : bin2beg
 Usage   : my $start = BGZFast::hts_utility::bin2beg($bin);
 Function: Get the region start position corresponding to a bin
 Returns : Start
 Args    : Bin number

=cut
sub bin2beg {
    return ( ($_[0] - 4681)<<14 )if($_[0] >= 4681);
    return ( ($_[0] -  585)<<17 )if($_[0] >=  585);
    return ( ($_[0] -   73)<<20 )if($_[0] >=   73);
    return ( ($_[0] -    9)<<23 )if($_[0] >=    9);
    return ( ($_[0] -    1)<<26 )if($_[0] >=    1);
    return 0;
}
#------------------------------------------------------------------------
=head2 bin2end

 Title   : bin2end
 Usage   : my $end = BGZFast::hts_utility::bin2end($bin);
 Function: Get the region end position corresponding to a bin
 Returns : End
 Args    : Bin number

=cut
sub bin2end {
    return (( ($_[0] - 4681)<<14 ) +    16383)if($_[0] >= 4681); #beg + ((1<<14) - 1)
    return (( ($_[0] -  585)<<17 ) +   131071)if($_[0] >=  585); #beg + ((1<<17) - 1)
    return (( ($_[0] -   73)<<20 ) +  1048575)if($_[0] >=   73); #beg + ((1<<20) - 1)
    return (( ($_[0] -    9)<<23 ) +  8388607)if($_[0] >=    9); #beg + ((1<<23) - 1)
    return (( ($_[0] -    1)<<26 ) + 67108863)if($_[0] >=    1); #beg + ((1<<26) - 1)
    return 536870911; #beg + ((1<<29) - 1)
}
#------------------------------------------------------------------------
=head2 bin2reg

 Title   : bin2reg
 Usage   : my ($start, $end) = BGZFast::hts_utility::bin2reg($bin);
 Function: Get the region [start, end) corresponding to a bin
 Returns : Start and end
 Args    : Bin number

=cut
sub bin2reg {
    if($_[0] >= 4681){
	$_[0] = ($_[0] - 4681)<<14;
	return ($_[0], $_[0] + 16383);
    }
    elsif($_[0] >=  585){
	$_[0] = ($_[0] - 585)<<17;
	return ($_[0], $_[0] + 131071);
    }
    elsif($_[0] >=   73){
	$_[0] = ($_[0] - 73)<<20;
	return ($_[0], $_[0] + 1048575);
    }
    elsif($_[0] >=    9){
	$_[0] = ($_[0] - 9)<<23;
	return ($_[0], $_[0] + 8388607);
    }
    elsif($_[0] >=    1){
	$_[0] = ($_[0] - 1)<<26;
	return ($_[0], $_[0] + 67108863);
    }
    else{
	return (0, 536870911);
    }
}
#------------------------------------------------------------------------
=head2 bin2parent

 Title   : bin2parent
 Usage   : my $parent = BGZFast::hts_utility::bin2parent($bin);
 Function: Identify immediate parent of a bin
 Returns : Bin numbers of parent
 Args    : Bin number

=cut
sub bin2parent {
    return (585 + (($_[0]-4681)>>3)) % 65536 if($_[0] >= 4681);
    return ( 73 + (($_[0]- 585)>>3)) % 65536 if($_[0] >=  585);
    return (  9 + (($_[0]-  73)>>3)) % 65536 if($_[0] >=   73);
    return (  1 + (($_[0]-   9)>>3)) % 65536 if($_[0] >=    9);
    return 0 if($_[0] >= 1);
    return undef;
}
#------------------------------------------------------------------------
=head2 bin2parents

 Title   : bin2parents
 Usage   : my @parents = BGZFast::hts_utility::bin2parents($bin);
 Function: Identify all parents of a bin
 Returns : List of bin numbers of parents
 Args    : Bin number

=cut
sub bin2parents {
    if($_[0] >= 4681){
        $_[0] -= 4681;
        return (0,
                (  1 + ($_[0]>>12)) % 65536,
                (  9 + ($_[0]>> 9)) % 65536,
                ( 73 + ($_[0]>> 6)) % 65536,
                (585 + ($_[0]>> 3)) % 65536);
    }
    elsif($_[0] >= 585){
        $_[0] -= 585;
        return (0,
                ( 1 + ($_[0]>>9)) % 65536,
                ( 9 + ($_[0]>>6)) % 65536,
                (73 + ($_[0]>>3)) % 65536);
    }
    elsif($_[0] >= 73){
        $_[0] -= 73;
        return (0,
                (1 + ($_[0]>>6)) % 65536,
                (9 + ($_[0]>>3)) % 65536);
    }
    elsif($_[0] >= 9){
        $_[0] -= 9;
        return (0,
                (1 + ($_[0]>>3)) % 65536);
    }
    elsif($_[0] >= 1){
        return 0;
    }
    else{
        return undef;
    }
}
#------------------------------------------------------------------------
=head2 bin2children

 Title   : bin2children
 Usage   : my @children = BGZFast::hts_utility::bin2children($bin);
 Function: Identify all children of a bin
 Returns : List of bin numbers of children
 Args    : Bin number

=cut
sub bin2children {
    my ($beg, $end) = bin2reg($_[0]);
    my $i = 0;
    my @list;
    if($_[0] < 1){
	push(@list, 1..8);
	$i+=8;
    }
    if($_[0] < 9){
	for(my $k = 9 + ($beg>>23); $k <= 9 + ($end>>23); $k++){$list[$i++] = $k}
    }
    if($_[0] < 73){
	for(my $k = 73 + ($beg>>20); $k <= 73 + ($end>>20); $k++){$list[$i++] = $k}
    }
    if($_[0] < 585){
	for(my $k = 585 + ($beg>>17); $k <= 585 + ($end>>17); $k++){$list[$i++] = $k}
    }
    if($_[0] < 4681){
	for(my $k = 4681 + ($beg>>14); $k <= 4681 + ($end>>14); $k++){$list[$i++] = $k}
    }

    return \@list;
}
#------------------------------------------------------------------------
#----                        FEATURE RELATED                         ----
#------------------------------------------------------------------------
=head2 feature_stat

 Title   : feature_stat
 Usage   : my ($length, $chr, $pos, $end, $next_chr, $next_pos) =
             BGZFast::hts_utility::feature_stat($string, $offset, 'VCF');

           my @stat = BGZFast::hts_utility::feature_stat($string,
                                                         $offset,
                                                         'GENERIC',
                                                         [1, 4, 5]);
 Function: Get basic attributes of a feature line
 Returns : Feature length, ref_name, beg, end, next_name, next_pos
 Args    : String, offset, index_type, [index columns in array ref]

 File Format   index_type INT   index_type TEXT
 -----------   --------------   ---------------
 Sam format    1                'SAM'
 Vcf format    2                'VCF'
 Generic       0                'GENERIC'
 Zero-based    65536            'ZERO'

For multipart features or features with breakpoints, the position of the
next feature/breakpoint are given in next_name and next_pos values. All
positions returned are zero based. The index_type value can either be
a text or an integer value. An index column array reference is only
required for generic and zero-based file types (see Tabix documentation
for more info).

=cut
sub feature_stat {
    my $ref  = (ref($_[0])) ? $_[0] : \ ($_[0]);
    my $boff = $_[1] || 0;
    my $type = $_[2] || 0;  #index type
    my $c    = $_[3];       #index cols
    my $max  = $_[4] || 2147483647; #index col max

    #find length and end offset
    my $eoff = index($$ref, "\n", $boff)+1;
    return undef if(!$eoff); #not a full feature
    my $length = $eoff-$boff;
    $eoff--; #set end offset before endline

    #convert text type to numerical type (condensed format)
    if(ord(substr($type, 0, 1)) > 57){
        $type = uc($type);
        if   ($type eq     'SAM'){$type =     1; $max = 6}
        elsif($type eq     'VCF'){$type =     2; $max = 8}
        elsif($type eq 'GENERIC'){$type =     0; die "ERROR: Missing index cols\n" if(!$c)}
        elsif($type eq    'ZERO'){$type = 65536; die "ERROR: Missing index cols\n" if(!$c)}
        else{die "ERROR: Invalid index type given: $type\n"}
    }
    elsif($type == 1){$max = 6} #fix max
    elsif($type == 2){$max = 8} #fix max
    elsif(!$c){die "ERROR: Missing index cols\n"}

    #get tabs to avoid using split
    my @t = ($boff); #tab positions (offset, length)
    for(my $i = 1; $i <= $max; $i++){
        $t[$i*2] = index($$ref, "\t", $t[$i-1])+1; #set offset
        $t[$i*2-1] = $t[$i*2] - $t[$i*2-2]-1; #set prev length
        if(!$t[$i*2] || $t[$i*2] > $eoff){ #fix if past end
            shift(@t); #drop partial
            $t[$i*2-1] = $eoff - $t[$i*2-2];
            last;
        }
        elsif($i == $max){ #drop partial
            shift(@t);
        }
    }

    #get stats for type
    my ($chr, $pos, $end, $next_chr, $next_pos);
    if($type == 1){ #SAM
        $chr = substr($$ref, $t[2*2], $t[2*2+1]); #col 3

        #get position
        if($chr eq '*'){ #unmapped
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
        $next_chr = substr($$ref, $t[6*2], $t[6*2+1]); #col 7
        $next_chr = $chr if($next_chr eq '=');
        $next_pos = ($next_chr eq '*') ? -1 : substr($$ref, $t[7*2], $t[7*2-1])-1; #col 8
    }
    elsif($type == 2){ #VCF
        $chr = substr($$ref, $t[0*2], $t[0*2+1]); #col 1
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
            $next_chr = $1;
            $next_pos = $2-1; #make zero based
            if(index($alt, ',', 0)+1){ #too complex
                warn "WARNING: Feature is too complex:\n".substr($$ref, $boff, $eoff-$boff)."\n";
                $next_pos   = -2;
                $next_chr = \0;
            }
        }
    }
    elsif($type == 0){ #Generic
        $chr = substr($$ref, $t[$c->[0]*2-2], $t[$c->[0]*2-1]);
        if($chr eq '*'){
            $pos = -1;
            $end = -1;
        }
        else{
            $pos = substr($$ref, $t[$c->[1]*2-2], $t[$c->[1]*2-1])-1; #-1 to make zero based
            $end = ($c->[2] && $c->[2] != $c->[1]) ?
                substr($$ref, $t[$c->[2]*2-2], $t[$c->[2]*2-1])-1 : $pos+1;
        }
        $next_pos = -1;
        $next_chr = '*';
    }
    elsif($type == 65536){ #Zero-based
        my $chr = substr($$ref, $t[$c->[0]*2-2], $t[$c->[0]*2-1]);
        if($chr eq '*'){
            $pos = -1;
            $end = -1;
        }
        else{
            $pos = substr($$ref, $t[$c->[1]*2-2], $t[$c->[1]*2-1]);
            $end = ($c->[2] && $c->[2] != $c->[1]) ?
                substr($$ref, $t[$c->[2]*2-2], $t[$c->[2]*2-1])-1 : $pos+1;
        }
        $next_pos = -1;
        $next_chr = '*';
    }
    else{
        die "ERROR: invalid index type given: $type\n"
    }

    return ($length, $chr, $pos, $end, $next_chr, $next_pos);
}
#------------------------------------------------------------------------
1;
