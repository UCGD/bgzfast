#!/usr/bin/perl

use forks;
use forks::shared;

use FindBin;
use lib "$FindBin::RealBin/../lib";

use strict;
use warnings;
use Carp;
use Perl::Unsafe::Signals;
use Getopt::Long qw(:config no_ignore_case);
use POSIX qw(ceil);
use Fcntl qw(:seek);
use BGZFast::vcf_manipulator;
use BGZFast::tbx_extra qw(tbx_index_vcf_chunk);

BEGIN{
    binmode(STDIN);
    binmode(STDOUT);
    select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
    select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately

    $SIG{INT}  = sub {exit(130)};
    $SIG{USR2} = sub {print STDERR Carp::longmess;}; #set signal handler

    mkdir("$ENV{HOME}/.Inline") if(! -d "$ENV{HOME}/.Inline");
}

#set globals
our $DEBUG;
our $VERSION = 0.0.1;

my $usage = "
Version: $VERSION
    Usage:   vcf_indexer.pl [options] <file>

    This script performs parallel tabix based indexing of vcf files.

Standard options:
   -C, --cpus        INT     Number of CPUs to use for indexing [1]
   -f, --force               overwrite existing index without asking

";

my $preset;
my $cpus = 1;
my $force;
GetOptions("preset|p=s"   => \$preset,
	   "cpus|C=i"     => \$cpus,
	   "force|f"      => \$force,
	   "DEBUG|debug"  => \$DEBUG,
	   "help|?"       => sub {print $usage; exit(0);}
    ) or die("Error in command line arguments\n");

my $file = shift;

#usage and initial error handling
if(!$file){
    print $usage;
    die "ERROR: Failure to specify file\n";
}
die "Error: Input file does not exists: $file\n" if(! -f $file);

#handle presets

#create tabix object
my $vcf = BGZFast::vcf_manipulator->new($file);
					
#check index file
if(-f $vcf->tbi_file && !$force){
    die "ERROR: Index file already exists, use --force\n";
}

#get window size for job splitting
my $size = (stat($file))[7];
my $window = ceil($size/ceil($size/33554432));
my $count = ceil($size/$window);
$cpus = $count if(ceil($count < $cpus));

#distribute data to threads
my $tbi;
$vcf->close_handle(); #prep for forking
my @todo :shared = (1..$count);

if($DEBUG){
    @todo =(1..$count);
    while(my $id = shift(@todo)){
	my $start_time = time();
	my $tbi2 = index_part($vcf, [$id]);
	$tbi = ($tbi) ? merge_index($tbi, $tbi2) : $tbi2;
	my $run_time = time()-$start_time;
	print STDERR "##RUNTIME $id: $run_time\n";
    }
}
else {
    for(my $i = 1; $i < $cpus; $i++){
	threads->create({'context' => 'scalar'}, \&index_part, $vcf, \@todo);
    }
    $tbi = index_part($vcf, \@todo);

    #gather index parts
    while(threads->list()){
	my @thr = threads->list(threads::joinable);
	if(!@thr){
	    sleep 0.01;
	    next;
	}

	foreach my $thr (@thr){
	    my $tbi2 = $thr->join;
	    $tbi = merge_index($tbi, $tbi2);
	}
    }
}

#fix missing bins in linear index
my $index = $tbi->indices;
for(my $i = 0; $i < @$index; $i++){
    my $lindex = $index->[$i]{LINEAR};
    for(my $l = $#$lindex; $l >= 0; $l--){
	$lindex->[$l] = [@{$lindex->[$l+1]}] if(!$lindex->[$l]);
    }
}

#collapse small bins onto parent if bin is too small
for(my $i = 0; $i < @$index; $i++){
    my $bindex = $index->[$i]{BINS};
    for(my $bin = $#$bindex; $bin > 0; $bin--){
	next if(!$bindex->[$bin]);

	#move chunks to parent
	my $p = BGZFast::vcf_manipulator->bin2parent($bin);
	if($bindex->[$bin][-1][2] - $bindex->[$bin][0][0] < 65536 && $bindex->[$p]){ #at least 64kb
	    my @chunks = (@{$bindex->[$p]}, @{$bindex->[$bin]});
	    @chunks = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @chunks;
	    $bindex->[$p] = \@chunks;
	    delete($bindex->[$bin]); #removes all trailing non-existing elements
	    $bin = @$bindex if($bin > @$bindex); #fix for potential truncation
	    next;
	}
	
	#merge chunks on same block
	my $chunks = $bindex->[$bin];
	for(my $j = 0; $j < $#$chunks; $j++){
	    my $k = $j+1;
	    next unless($chunks->[$k][0] <= $chunks->[$j][2]);
	    if($chunks->[$k][2] > $chunks->[$j][2] ||
	       ($chunks->[$k][2] == $chunks->[$j][2] && $chunks->[$k][3] > $chunks->[$j][3])
		){
		$chunks->[$j][2] = $chunks->[$k][2];
		$chunks->[$j][3] = $chunks->[$k][3];
	    }
	    splice(@$chunks,$k,1); #remove $k
	    $j--; #step back to account for splice
	}
    }
}

#write the index
$vcf->{TBI} = $tbi; #override any internal index
$vcf->write_tbi();

#-------------------------------------------------------------------------------
#--------------------------------- SUBROUTINES ---------------------------------
#-------------------------------------------------------------------------------

#convenience method to see if I am in a thread or not
{my $tid; BEGIN{$tid = threads->tid if exists $INC{'threads.pm'}};
sub is_thread {
    my $is_thread = 1 if(defined($tid) && $tid != threads->tid);
    return $is_thread;
}}

sub index_part {
    my $vcf = shift;
    my $todo = shift;

    #fix thread signalling
    if(is_thread){
        select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately
	select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
	
	$SIG{'__DIE__'} = sub {print STDERR  "$_[0]\n"; exit(255)};
        $SIG{INT} = sub {exit(130)};
    }
    
    #make empty index hash
    my $tbi = $vcf->__tbi_struct;
    $tbi->type = $vcf->{IDX_TYPE}; # 2 for vcf
    $tbi->col  = [@{$vcf->{IDX_COL}}]; # 1,2, 0 for vcf
    $tbi->meta = $vcf->{IDX_META}; # '#' for vcf
    $tbi->skip = $vcf->{IDX_SKIP}; # 0 for vcf

    #read in data to index
    my $f_size = $vcf->__file->length;
    while(my $id = shift @$todo){
	my $w_beg = ($id-1)*$window;
	my $w_end = ($id-0)*$window;
	$w_end = $f_size if($w_end > $f_size);

	my ($head_foff, $head_voff) = $vcf->smart_seek($w_beg);
	my ($tail_foff, $tail_voff) = $vcf->smart_seek($w_end);

	my $tbi_string = tbx_index_vcf_chunk($vcf->file, $head_foff, $head_voff, $tail_foff, $tail_voff);
	my $tbi2 = $vcf->__string2tbi($tbi_string);
	$tbi = merge_index($tbi, $tbi2);
    }
    
    return $tbi;
}

sub add_to_index {
    my $tbi = shift;
    my $dref = shift; #uncompressed
    my $offsets = shift; #(buf_off, length, file_off, vblock_off)

    #get format information
    my $index  = $tbi->indices;
    my $names  = $tbi->names;
    my $name2tbid = $tbi->name2tbid;
    my $max = 8; #max columns we need to read is 8 for vcf files

    #process each feature
    my @tpos;
    my $lbeg = 0; #line begin offset
    my $lend = index($$dref, "\n", $lbeg)+1; #line end offset
    while($lend){
	my ($seq, $B, $E) = @{feature_pos($$dref, $lbeg, $lend)};

	#get real and virtual offsets
	shift(@$offsets) while($lbeg >= $offsets->[0][0] + $offsets->[0][1]);
	my ($foffB, $voffB) = ($offsets->[0][2], $offsets->[0][3] + $lbeg-$offsets->[0][0]);
	shift(@$offsets) while($lend > $offsets->[0][0] + $offsets->[0][1]);
	my ($foffE, $voffE) = ($offsets->[0][2], $offsets->[0][3] + $lend-$offsets->[0][0]);
	($foffE, $voffE) = ($offsets->[0][2]+$offsets->[0][4], 0)
	    if($lend == $offsets->[0][0] + $offsets->[0][1]); #pos at end of virtual offsets

	#add to indexes
	if($seq eq '*'){ #unmapped only get counts
	    $tbi->n_no_coor++; #itterate count
	}
	else{ #mapped reads go into index
	    #initialize reference
	    my $tid = $name2tbid->{$seq};
	    if(!defined($tid)){
		push(@$names, $seq);
		$tid = $name2tbid->{$seq} = $#$names;
		$index->[$tid] = {BINS => [], LINEAR => [], PSEUDO => []};
	    }
	    
	    #get index parts
	    my $bindex = $index->[$tid]{BINS};
	    my $lindex = $index->[$tid]{LINEAR};
	    my $pseudo = $index->[$tid]{PSEUDO};

	    #update linear index
	    my $bin = BGZFast::vcf_manipulator->reg2bin($B, $E);
	    if($bin >= 4681){
		$lindex->[$bin-4681] = [$foffB, $voffB] if(!$lindex->[$bin-4681]);
	    }
	    else{
		my $lins = BGZFast::vcf_manipulator->reg2lins($B, $E);
		foreach my $l (@$lins){
		    $lindex->[$l] = [$foffB, $voffB] if(!$lindex->[$l]);
		}
	    }

	    #update bin index
	    if(!$bindex->[$bin]){ #initialize
		$bindex->[$bin] = [[$foffB, $voffB, $foffE, $voffE]];
	    }
	    elsif($foffB <= $bindex->[$bin][-1][2]){ #ends near same block as current
		$bindex->[$bin][-1][2] = $foffE; #extend physical end
		$bindex->[$bin][-1][3] = $voffE; #extend virtual end
	    }
	    else{ #make new chunk
		push(@{$bindex->[$bin]}, [$foffB, $voffB, $foffE, $voffE]);
	    }

	    #update pseudobin
	    if(!@$pseudo){ #initialize
		$pseudo->[0] = [$foffB, $voffB, $foffE, $voffE];
		$pseudo->[1] = [0, 0]; #mapped_count, placed_unmapped_count
	    }
	    else{
		$pseudo->[0][2] = $foffE; #extend physical end
		$pseudo->[0][3] = $voffE; #extend virtual end
	    }

	    #add counts to pseudobin
	    $pseudo->[1][0]++; #itterate mapped count
	}
	    
	#reset for next feature
	$lbeg = $lend;
	$lend = index($$dref, "\n", $lbeg)+1;
    }

    #trim already processed
    if($lbeg){
	substr($$dref, 0, $lbeg, '');
	if(@$offsets = grep {$lbeg < $_->[0]+$_->[1]} @$offsets){
	    $_->[0] -= $lbeg foreach(@$offsets);
	    $offsets->[0][1] += $offsets->[0][0];
	    $offsets->[0][3] -= $offsets->[0][0];
	    $offsets->[0][0] = 0;
	}
    }

    return $lbeg; #start position of last partial feature
}

sub merge_index {
    my $tbi1 = shift;
    my $tbi2 = shift;

    return $tbi1 if(!$tbi2);
    return $tbi2 if(!$tbi1);

    #get proper order
    my %ps;
    my $names1 = $tbi1->names;
    my $names2 = $tbi2->names;
    for(my $i = 0; $i < @$names1; $i++){
	$ps{$names1->[$i]} = $tbi1->indices($i)->{PSEUDO}[0];
    }
    for(my $i = 0; $i < @$names2; $i++){
	$ps{$names2->[$i]} = $tbi2->indices($i)->{PSEUDO}[0];
    }
    my @order = sort {$ps{$a}[0] <=> $ps{$b}[0] || $ps{$a}[1] <=> $ps{$b}[1]} keys(%ps);

    #make empty index hash
    my $tbi = BGZFast::vcf_manipulator->__tbi_struct;
    $tbi->type = $tbi1->type;
    $tbi->col  = [@{$tbi1->col}];
    $tbi->meta = $tbi1->meta;
    $tbi->skip = $tbi1->skip;
    $tbi->names = \@order;
    $tbi->name2tbid = {map {$order[$_] => $_} (0..$#order)};
    $tbi->n_no_coor = $tbi1->n_no_coor + $tbi2->n_no_coor;

    #iterate through references and add to index
    for(my $i = 0; $i < @order; $i++){
	#get ids from other indexes
	my $seq = $order[$i];
	my $id1 = $tbi1->name2tbid($seq);
	my $id2 = $tbi2->name2tbid($seq);

	#does not exist
	if(!defined($id1)){
	    $tbi->indices->[$i] = $tbi2->indices($id2);
	    next;
	}

	if(!defined($id2)){
	    $tbi->indices->[$i] = $tbi1->indices($id1);
	    next;
	}

	#initialize
	$tbi->indices->[$i] = $tbi1->indices($id1);

	#merge linear index
	my $lindex1 = $tbi1->indices($id1)->{LINEAR};
	my $lindex2 = $tbi2->indices($id2)->{LINEAR};
	for(my $l = 0; $l < @$lindex2; $l++){
	    if(!$lindex2->[$l]){
		next;
	    }
	    elsif(!$lindex1->[$l]){
		$lindex1->[$l] = $lindex2->[$l];
		next;
	    }
	    elsif($lindex2->[$l][0] > $lindex1->[$l][0]){
		next;
	    }
	    elsif($lindex2->[$l][0] < $lindex1->[$l][0]){
		$lindex1->[$l] = $lindex2->[$l];
		next;
	    }
	    elsif($lindex2->[$l][1] < $lindex1->[$l][1]){ #implicit ($lindex2->[$l][0] == $lindex1->[$l][0])
		$lindex1->[$l] = $lindex2->[$l];
		next;
	    }
	}
	
	#merge bin index
	my $bindex1 = $tbi1->indices($id1)->{BINS};
	my $bindex2 = $tbi2->indices($id2)->{BINS};
	for(my $bin = 0; $bin < @$bindex2; $bin++){
	    if(!$bindex2->[$bin]){ #nothing
		next;
	    }
	    elsif(!$bindex1->[$bin]){ #replace
		$bindex1->[$bin] = $bindex2->[$bin];
		next;
	    }
	    else{ #merge
		my @chunks = (@{$bindex1->[$bin]}, @{$bindex2->[$bin]});
		@chunks = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @chunks;
		for(my $j = 0; $j < $#chunks; $j++){
		    my $k = $j+1;
		    next unless($chunks[$k][0] <= $chunks[$j][2]);
		    if($chunks[$k][2] > $chunks[$j][2] ||
		       ($chunks[$k][2] == $chunks[$j][2] && $chunks[$k][3] > $chunks[$j][3])
			){
			$chunks[$j][2] = $chunks[$k][2];
			$chunks[$j][3] = $chunks[$k][3];
		    }
		    splice(@chunks,$k,1); #remove $k
		    $j--; #step back to account for splice
		}
		$bindex1->[$bin] = \@chunks;
		next;
	    }
	}
	
	#merge pseudobin
        my $pseudo  = $tbi->indices($i)->{PSEUDO};
        my $pseudo2 = $tbi2->indices($id2)->{PSEUDO};
	if($pseudo2->[0][0] < $pseudo->[0][0] ||
	   ($pseudo2->[0][0] == $pseudo->[0][0] && $pseudo2->[0][1] < $pseudo->[0][1])
	){ #grow
	    $pseudo->[0][0] = $pseudo2->[0][0];
	    $pseudo->[0][1] = $pseudo2->[0][1];
	}
	if($pseudo2->[0][2] > $pseudo->[0][2] ||
	   ($pseudo2->[0][2] == $pseudo->[0][2] && $pseudo2->[0][3] > $pseudo->[0][3])
	){ #grow
	    $pseudo->[0][2] = $pseudo2->[0][2];
	    $pseudo->[0][3] = $pseudo2->[0][3];
	}
	$pseudo->[1][0] += $pseudo2->[1][0]; #iterate
	$pseudo->[1][1] += $pseudo2->[1][1]; #iterate
    }

    return $tbi;
}
