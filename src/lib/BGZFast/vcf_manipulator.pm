#------------------------------------------------------------------------
#----                   BGZFast::vcf_manipulator                     ----
#------------------------------------------------------------------------
package BGZFast::vcf_manipulator;

use strict;
use warnings;
use Compress::Zlib;
use Fcntl qw(:DEFAULT :seek); #SEEK_SET, SEEK_CUR, SEEK_END

use parent qw(BGZFast::hts_manipulator);

our $VERSION = '0.01';

#------------------------------------------------------------------------
#--------------       OBJECT INITIALIZATION METHODS        --------------
#------------------------------------------------------------------------

sub new {
    return shift->SUPER::new(@_);
}

sub __init {
    my $self = shift;
    my @args = @_;

    #initialize base class parameters (returns remaining params)
    my $param = $self->SUPER::__init(@args);

    #force index parameters
    $self->{IDX_TYPE} = 2; #VCF code
    $self->{IDX_COL}  = [1,2,0];
    $self->{IDX_META} = '#';
    $self->{IDX_SKIP} = 0;

    #parse header
    $self->__parse_header();

    return $param;
}

#------------------------------------------------------------------------
#--------------    OBJECT SPECIFIC METHODS TO OVERRIDE     --------------
#------------------------------------------------------------------------

#parses info from feature text (should be overridden by child module)
sub __fill_feature_ref {
    my $self = shift;
    my $feat = shift;

    my @F = split(/\t/, $feat->data);
    chomp($F[-1]);

    $feat->refid = $self->name2seqid($F[0]);
    $feat->pos = $F[1]-1; #make zero based
    $feat->end = $F[1]+length($F[3])-1; #make zero based
    $feat->next_ref_id = -1;
    $feat->next_pos = -1;
    
    #add breakpoints for SVs
    if($F[4] =~ /[\[\]]\s*([^\[\]\:\s]+)\s*\:\s*([^\[\]\:\s]+)\s*[\[\]]/){
	$feat->next_pos = $2-1; #make zero based
	$feat->next_ref_id = $self->name2seqid($1);
    }

    return;
}

#parses header reference text (should be overridden by child module)
sub __fill_header_ref {
    my $self = shift;
    my $header = shift;

    my $ref = $header->ref;
    my $name2seqid = $header->name2seqid;
    my $name2length = $header->name2length;
    
    #parse reference lines
    my $i = 0;
    foreach my $line (split(/\n/, $header->text)){
	next unless($line =~ /^##contig=/);
	my ($name) = $line =~ /[\,\<]ID=([^\,\>]+)/;
	my ($l_ref) = $line =~ /[\,\<]length=([^\,\>]+)/;
	
	$ref->[$i]{NAME} = $name;
	$ref->[$i]{LENGTH} = $l_ref;
	$name2seqid->{$name} = $i;
	$name2length->{$name} = $l_ref;
	$i++;
    }
    $header->n_ref = $i;

    return;
}

#------------------------------------------------------------------------
#--------------              GENERIC METHODS               --------------
#------------------------------------------------------------------------

sub seqid2name {
    my $self = shift;
    my $id = shift;
    my $header = $self->__header;

    if(! defined($header->ref($id))){
	if($self->tbi){
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

sub seqid2seqlength {
    my $self = shift;
    my $id = shift;
    my $header = $self->__header;

    if(! defined($header->ref($id))){
	warn "WARNING: No contig exists with seqID $id\n";
	return undef;
    }

    return $header->ref($id)->{LENGTH};
}

sub name2seqid {
    my $self = shift;
    my $name = shift;
    my $header = $self->__header;

    if(!defined($header->name2seqid($name))){
	if($self->tbi && !$self->tbi->{__HEADERFILL}){
	    my $tbi = $self->tbi;
	    my $names = $tbi->names;
	    for(my $i = 0; $i < @$names; $i++){
		$header->name2seqid($names->[$i]) ||= $i;
	    }
	    $self->tbi->{__HEADERFILL} = 1; #flag to only use IDX to fill header once

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

sub name2seqlength {
    my $self = shift;
    my $name = shift;
    my $header = $self->__header;

    if(! defined($header->name2length($name))){
	warn "WARNING: No contig exists with name $name\n";
	return undef;
    }

    return $header->name2length($name);
}

1;
