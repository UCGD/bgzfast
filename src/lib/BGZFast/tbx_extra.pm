#!/usr/bin/perl

package BGZFast::tbx_extra;

use base qw(Exporter);
@EXPORT_OK = qw(tbx_index_vcf_chunk);

use strict;
use warnings;
use Cwd;
use File::Which;
use Inline;
use vars qw($CODE $INC $HTS_INC $HTS_LIB);

our $VERSION = '0.01';

#C code is compiled and included here
if($HTS_INC && $HTS_LIB){
    Inline->bind(C => $CODE,
		 clean_after_build => 0,
		 NAME => 'BGZFast::tbx_extra',
		 VERSION => ($INC{'BGZFast/ConfigData.pm'}) ? $VERSION : undef,
		 LIBS => "-L$HTS_LIB -lhts",
		 INC => "-I$HTS_INC -I$INC");
}
else{
    warn "WARNING: ".__PACKAGE__." not configured correctly. Can't find hts_lib\n";
}

sub tbx_index_vcf_chunk { return __tbx_index_vcf_chunk(@_); }

BEGIN {
    eval "require BGZFast::ConfigData";

    if($INC{'BGZFast/ConfigData.pm'}){
	$HTS_INC = BGZFast::ConfigData->config('HTS_INC');
	$HTS_LIB = BGZFast::ConfigData->config('HTS_LIB');
    }
    elsif(my $tabix = File::Which::which('tabix')){
	(my $HTS = $tabix) =~ s/(?:bin\/)?tabix$//;
	$HTS_INC = "$HTS/include";
	$HTS_LIB = "$HTS/lib";
    }

    ($INC = Cwd::abs_path(__FILE__)) =~ s/tbx_extra.pm$//;
    
    $CODE = <<END;
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "tbx_extra.h"

SV* __tbi2string(tbx_t *tbx) {
    hts_idx_t* idx = tbx->idx;
    size_t max = 100000;
    char* string = (char*) malloc(100000*sizeof(char));
    size_t len = 0;

    memcpy(&string[len], "TBI\\1", 4);
    len += 4;

    int32_t i, j;
    int nids = idx->n;
    for (i = nids = 0; i < idx->n; ++i) {
	if (idx->bidx[i]) nids++;
    }

    memcpy(&string[len], &nids, 4);
    len += 4;

    if(max < len+idx->l_meta){max += 100000; string = (char*) realloc(string, max);}
    memcpy(&string[len], idx->meta, idx->l_meta);
    len += idx->l_meta;

    for (i = 0; i < idx->n; ++i) {
	khint_t k;
	bidx_t *bidx = idx->bidx[i];
	lidx_t *lidx = &idx->lidx[i];

	// write binning index
	if (nids == idx->n || bidx){
            int32_t n_bin = bidx ? kh_size(bidx) : 0;

            if(max < len+4){max += 100000; string = (char*) realloc(string, max);}
            memcpy(&string[len], &n_bin, 4);
            len += 4;
	}

        if (bidx){
            for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
		if (kh_exist(bidx, k)) {
                    bins_t *p = &kh_value(bidx, k);

                    if(max < len+4){max += 100000; string = (char*) realloc(string, max);}
                    memcpy(&string[len], &kh_key(bidx, k), 4);
                    len += 4;

                    if(max < len+4){max += 100000; string = (char*) realloc(string, max);}
                    memcpy(&string[len], &p->n, 4);
                    len += 4;

                    for (j = 0; j < p->n; ++j) {
			if(max < len+8){max += 100000; string = (char*) realloc(string, max);}
			memcpy(&string[len], &p->list[j].u, 8);
			len += 8;

			if(max < len+8){max += 100000; string = (char*) realloc(string, max);}
			memcpy(&string[len], &p->list[j].v, 8);
			len += 8;
                    }
		}
            }
	}

        // write linear index
	if(max < len+4){max += 100000; string = (char*) realloc(string, max);}
	memcpy(&string[len], &lidx->n, 4);
	len += 4;

	for (j = 0; j < lidx->n; ++j){
            if(max < len+8){max += 100000; string = (char*) realloc(string, max);}
            memcpy(&string[len], &lidx->offset[j], 8);
            len += 8;
	}
    }

    if(max < len+8){max += 100000; string = (char*) realloc(string, max);}
    memcpy(&string[len], &idx->n_no_coor, 8);
    len += 8;

    SV* ret = newSVpvn(string, len);
    free(string);

    return ret;
}

SV* __tbx_index_vcf_chunk(SV* name, unsigned long f_beg, unsigned long v_beg, unsigned long f_end, unsigned long v_end) {
    char* fn = SvPV_nolen(name);
    uint64_t c_beg = (f_beg << 16)|(v_beg & 0xFFFF);
    uint64_t c_end = (f_end << 16)|(v_end & 0xFFFF);

    BGZF *fp;
    fp = bgzf_open(fn, "r");
    bgzf_seek(fp, c_beg, SEEK_SET);

    //////////make index here
    tbx_t *tbx;
    kstring_t str;
    int stat, first = 0, min_shift = 14, n_lvls = 5, fmt = HTS_FMT_TBI;
    int64_t lineno = 0;
    uint64_t cur_off;
    uint64_t last_off = bgzf_tell(fp);
    tbx_intv_t intv;
    int64_t max_ref_len = 0;

    str.s = 0; str.l = str.m = 0;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx) goto fail;
    tbx->conf = tbx_conf_vcf;
    while ((stat = bgzf_getline(fp, '\\n', &str)) >= 0) {
	++lineno;
	cur_off = bgzf_tell(fp);

	if (str.s[0] == tbx->conf.meta_char) {
            last_off = cur_off;
            continue;
	}
        if (first == 0) {
            tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);
            if (!tbx->idx) goto fail;
            first = 1;
	}
        stat = get_intv(tbx, &str, &intv, 1);
	if (stat < -1) goto fail;  // Out of memory
	if (stat < 0) { if(c_end <= cur_off) break; continue; } // Skip unparsable lines
	if (hts_idx_push(tbx->idx, intv.tid, intv.beg, intv.end, cur_off, 1) < 0) goto fail;
	if(c_end <= cur_off) break; //chunk end
    }
    if (stat < -1) goto fail;
    if ( !tbx->idx ) tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);   // empty file
    if (!tbx->idx) goto fail;
    if ( !tbx->dict ) tbx->dict = kh_init(s2i);
    if (!tbx->dict) goto fail;
    if (hts_idx_finish_chunk(tbx->idx, bgzf_tell(fp)) != 0) goto fail;
    if (tbx_set_meta(tbx) != 0) goto fail;
    free(str.s);
    //////////

    bgzf_close(fp);

    SV* ret = __tbi2string(tbx);
    tbx_destroy(tbx);

    return ret;

 fail:
    croak("ERROR: tbx_index_vcf_chunk failed");
    return NULL;
}
	
END
}

1;
