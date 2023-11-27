#!/usr/bin/perl

#-------------------------------------------------------------------
package Bio::DB::FAI::Stream;
use base qw(Tie::Handle Bio::DB::SeqI);

sub new {
    my ($class, $db) = @_;
    my $key = $db->FIRSTKEY;
    return bless {
        db  => $db,
        key => $key
    }, $class;
}

sub next_seq {
    my $self = shift;
    my ($key, $db) = @{$self}{'key', 'db'};
    return if not defined $key;
    my $value = $db->get_Seq_by_id($key);
    $self->{key} = $db->NEXTKEY($key);
    return $value;
}

sub TIEHANDLE {
    my ($class, $db) = @_;
    return $class->new($db);
}

sub READLINE {
    my $self = shift;
    return $self->next_seq;
}

#-------------------------------------------------------------------
package Bio::DB::FAI;

use strict;
use warnings;
use Cwd;
use Inline;
use File::Spec;
use File::Temp qw(tempdir);
use base qw(Bio::DB::Fasta);
use vars qw($CODE $INC);

our $VERSION = '0.01';

#new method is implemented in C code below

#methods implemented in perl
sub seq {
    my ($self,$id,$start,$stop) = @_;

    if(!$start & !$stop && $id =~ /^(.+):([\d_]+)(?:,|-|\.\.)([\d_]+)$/ && !$self->exists($id)) {
	($id,$start,$stop) = ($1,$2,$3);
	$start =~ s/_//g;
	$stop =~ s/_//g;
	return $self->subseq($id,$start,$stop);
    }

    $start = 1 if($stop && !$start);
    $stop = $self->length($id) if($start && !$stop);

    my $reversed;
    if ($stop && $start > $stop) {
	($start,$stop) = ($stop,$start);
	$reversed++;
    }

    my $region = $id;
    $region   .= ":$start" if($start);
    $region   .= "-$stop"  if($stop);

    my $data =  $self->fetch($region);
    if ($reversed) {
	$data = reverse($data);
	$data =~ tr/gatcGATC/ctagCTAG/;
    }

    return $data;
}

sub subseq { return shift->seq(@_); }

sub ids { return shift->get_all_primary_ids(); }

sub get_Seq_by_id {
    my $self = shift;
    my $id   = shift;
    return unless $self->exists($id);
    return Bio::PrimarySeq::Fasta->new($self,$id);
}

sub alphabet {
    my $self = shift;
    my $id   = shift;
    
    local $_ = $self->seq($id);
    return ; /^[acgtrykmswbdhvnxACGTRYKMSWBDHVNX-]+$/  ? 'dna'
	   : /^[acgurykmswbdhvnxACGURYKMSWBDHVNX-]+$/  ? 'rna'
	   : /^[a-zA-Z*-]+$/  ? 'protein'
	   : '';
}

sub path { 
    my $self  = shift;
    return $self->filename;
}

sub file {
    my $self = shift;
    my $id = shift;

    return $self->filename if($self->exists($id));
}

sub index_name {
    my $self  = shift;
    my $path  = shift || $self->path;

    return "$path.fai";
}

sub get_PrimarySeq_stream {
    my $self = shift;
    return Bio::DB::FAI::Stream->new($self);
}

#inherited unsupported methods
sub newFh { die "ERROR: Bio::DB::FAI does not support initialization from filehandles\n"; }
sub dbmargs { die "ERROR: Bio::DB::FAI does not support dbmargs\n"; }
sub index_dir { die "ERROR: Bio::DB::FAI does not support indexing directories\n"; }
sub index_files { die "ERROR: Bio::DB::FAI does not support indexing multiple files\n"; }
sub offset { die "ERROR: Bio::DB::FAI does not support offset method\n"; }
sub header_offset { die "ERROR: Bio::DB::FAI does not support header_offset\n"; }
sub strlen { die "ERROR: Bio::DB::FAI does not support offset method\n"; }

#tie methods come here
sub TIEHASH { shift->new(@_); }
sub FETCH { shift->subseq(@_); }
sub STORE { shift->throw("Read-only database"); }
sub DELETE { shift->throw("Read-only database"); }
sub CLEAR { shift->throw("Read-only database"); }
sub EXISTS { shift->exists(@_); }
sub FIRSTKEY { shift->_firstkey(); }
sub NEXTKEY { shift->_nextkey(); }

#C code is compiled and included here
Inline->bind(C => $CODE,
	     clean_after_build => 0,
	     NAME => 'Bio::DB::FAI',
	     VERSION => $VERSION,
	     INC => "-I$INC");

BEGIN {
    ($INC = Cwd::abs_path(__FILE__)) =~ s/FAI.pm$//;
    $CODE = <<END;

#include <stdio.h>
#include "FAI.h"

SV* new(char* class, char* filename) {
    faidx_t* fai;
    SV*      obj_ref = newSViv(0);
    SV*      obj = newSVrv(obj_ref, class);

    fai = fai_load(filename);
    sv_setiv(obj, (IV)fai);
    SvREADONLY_on(obj);

    return obj_ref;
}

void index_file(SV* obj_ref, char* filename, int force) {
    faidx_t* old = (faidx_t*)SvIV(SvRV(obj_ref));
    fai_destroy(old);

    if(force != 0)
	unlink(filename);

    faidx_t* new;
    new = fai_load(filename);

    SV* obj = SvRV(obj_ref);
    SvREADONLY_off(obj);
    sv_setiv(obj, (IV)new);
    SvREADONLY_on(obj);

    return;
}

int exists(SV* obj, char* reg){
    faidx_t* fai = (faidx_t*)SvIV(SvRV(obj));
    char     *seq;
    return faidx_exists(fai,reg);
}

SV* fetch(SV* obj, char* reg){
    faidx_t* fai = (faidx_t*)SvIV(SvRV(obj));
    char     *seq;
    int       len;
    SV*       ret = NULL;

    seq = fai_fetch(fai,reg,&len);
    if(seq != NULL){
	ret = newSVpvn(seq,len);
	free((void*)seq);
    }

    return ret;
}

SV* header(SV* obj, char* reg){
    faidx_t* fai = (faidx_t*)SvIV(SvRV(obj));
    char     *head;
    SV*       ret = NULL;

    head = faidx_fetch_header(fai,reg);
    if(head != NULL){
	ret = newSVpv(head,0);
	free((void*)head);
    }

    return ret;
}

AV* get_all_primary_ids(SV* obj){
    faidx_t* fai;
    char**   keys;
    AV*      ret;

    fai  = (faidx_t*)SvIV(SvRV(obj));
    keys = faidx_all_keys(fai);
    ret = newAV();

    int i = 0;
    for(i=0; i<fai->n; i++){
	SV* id = NULL;
	id = newSVpv(keys[i],0);
	av_push(ret, id);
    }

    return ret;
}

char* filename(SV* obj){
    faidx_t* fai;
    fai  = (faidx_t*)SvIV(SvRV(obj));

    return fai->filename;
}

SV* length(SV* obj, char* reg){
    faidx_t* fai = (faidx_t*)SvIV(SvRV(obj));
    int64_t i64 = faidx_seq_length(fai,reg);
    SV* ret = NULL;

    if(i64 == -1) return ret;

    char str[21];
    int svlen = 0;
    while (i64) {
	str[svlen++] = (i64 % 10) + '0';
	i64 /= 10;
    }
    if (svlen) {
	ret = newSV(svlen);
        char *pv = SvPVX(ret);
        SvPOK_on(ret);
        SvCUR_set(ret, svlen);
	int i;
        for (i = svlen; i--;) *(pv++) = str[i];
    }
    else {
        ret = newSVpvs("0");
    }

    return ret;
}

char* _firstkey(SV* obj){
    faidx_t* fai;
    fai  = (faidx_t*)SvIV(SvRV(obj));
    char *key;

    fai->iter = 0;
    if(fai->iter < fai->n)
	key = fai->name[fai->iter++];
    
    return key;
}

char* _nextkey(SV* obj){
    faidx_t* fai;
    fai  = (faidx_t*)SvIV(SvRV(obj));
    char *key;

    if(fai->iter < fai->n)
	key = fai->name[fai->iter++];
    
    return key;
}

void DESTROY(SV* obj){
    faidx_t* fai = (faidx_t*)SvIV(SvRV(obj));
    fai_destroy(fai);
}
	
END
}

1;
