#!/usr/bin/perl

use FindBin;
use lib "$FindBin::RealBin/../lib";

use strict;
use warnings;
use Carp;
use Perl::Unsafe::Signals;
use BGZFast::vcf_manipulator;

my $vcf = BGZFast::vcf_manipulator->new('/scratch/general/pe-nfs1/u0045039/ptabix_test/NA12878.g.vcf.gz');
				
my $string1 = `gunzip -dc /scratch/general/pe-nfs1/u0045039/ptabix_test/NA12878.g.vcf.gz.sentieon`;
my $string2 = `gunzip -dc /scratch/general/pe-nfs1/u0045039/ptabix_test/NA12878.g.vcf.gz.tbi`;

#for(my$i=0; $i<length($string1); $i++){
#    my $c1 = substr($string1, $i, 1);
#    my $c2 = substr($string2, $i, 1);
#
#    print "string mismatch at $i\n" if($c1 ne $c2);
#}

my $tbi1 = $vcf->__string2tbi($string1);
exit();
my $tbi2 = $vcf->__string2tbi($string2);

my $string1b = $vcf->__tbi2string($tbi1);
my $string2b = $vcf->__tbi2string($tbi2);

for(my $i=0; $i < @{$tbi1->[7]}; $i++){
    my $bindex1 = $tbi1->[7][$i]{BINS};
    my $lindex1 = $tbi1->[7][$i]{LINEAR};
    my $pseudo1 = $tbi1->[7][$i]{PSEUDO};
    
    my $bindex2 = $tbi2->[7][$i]{BINS};
    my $lindex2 = $tbi2->[7][$i]{LINEAR};
    my $pseudo2 = $tbi2->[7][$i]{PSEUDO};

    #check pseudo
    my $ok = compare($pseudo1, $pseudo2);
    warn "pseudo mismatch at $i\n" if(!$ok);

    $ok = compare($bindex1, $bindex2);
    warn "bindex mismatch at $i\n" if(!$ok);

    $ok = compare($lindex1, $lindex2);
    warn "lindex mismatch at $i\n" if(!$ok);
    
    next;
}

exit(0);


sub compare {
    my $A = shift;
    my $B = shift;

    return 0 if(ref($A) ne ref($B));

    if(ref($A) eq ''){
	return 1 if(!defined($A) && !defined($B));
	return 0 if(!defined($A) || !defined($B));
	return ($A eq $B) ? 1 : 0;
    }
    elsif(ref($A) eq 'ARRAY'){
	return 0 if (@$A != @$B);

	for(my $i=0; $i<@$A; $i++){
	    return 0 if(!compare($A->[$i], $B->[$i]));
	}
	
	return 1;
    }
    elsif(ref($A) eq 'HASH'){
	my @keys_A = sort {$a cmp $b} keys(%$A);
	my @keys_B = sort {$a cmp $b} keys(%$B);
	
	return 0 if(!compare(\@keys_A, \@keys_B));
	
	foreach my $key (@keys_A){
	    return 0 if(!compare($A->{$key}, $B->{$key}));
	}

	return 1;
    }
    elsif(ref($A) eq 'SCALAR'){
	return compare($$A, $$B);
    }
    else{
	die "ERROR: Unsuported type:".ref($A);
    }

    return 1;
}

