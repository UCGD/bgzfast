
/// Format-neutral I/O, indexing, and iterator API functions.
/*
    Copyright (C) 2012-2021 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.
    Portions copyright (C) 2003-2006, 2008-2010 by Heng Li <lh3@live.co.uk>

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */


/// Tabix API functions.
/*
    Copyright (C) 2009, 2012-2015, 2019 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef TBX_EXTRA_H
#define TBX_EXTRA_H

#include "htslib/khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)
  
typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    hts_pos_t n, m;
    uint64_t *offset;
} lidx_t;

struct hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    int tbi_n, last_tbi_tid;
    struct {
        uint32_t last_bin, save_bin;
	hts_pos_t last_coor;
        int last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

static int tbx_set_meta(tbx_t *tbx)
{
    int i, l = 0, l_nm;
    uint32_t x[7];
    char **name;
    uint8_t *meta;
    khint_t k;
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;

    memcpy(x, &tbx->conf, 24);
    name = (char**)malloc(sizeof(char*) * kh_size(d));
    if (!name) return -1;
    for (k = kh_begin(d), l = 0; k != kh_end(d); ++k) {
        if (!kh_exist(d, k)) continue;
        name[kh_val(d, k)] = (char*)kh_key(d, k);
        l += strlen(kh_key(d, k)) + 1; // +1 to include '\0'
    }
    l_nm = x[6] = l;
    meta = (uint8_t*)malloc(l_nm + 28);
    if (!meta) { free(name); return -1; }
    if (ed_is_big())
        for (i = 0; i < 7; ++i)
            x[i] = ed_swap_4(x[i]);
    memcpy(meta, x, 28);
    for (l = 28, i = 0; i < (int)kh_size(d); ++i) {
        int x = strlen(name[i]) + 1;
        memcpy(meta + l, name[i], x);
        l += x;
    }
    free(name);
    hts_idx_set_meta(tbx->idx, l, meta, 0);
    return 0;
}

static inline int get_tid(tbx_t *tbx, const char *ss, int is_add)
{
    khint_t k;
    khash_t(s2i) *d;
    if (tbx->dict == 0) tbx->dict = kh_init(s2i);
    if (!tbx->dict) return -1; // Out of memory
    d = (khash_t(s2i)*)tbx->dict;
    if (is_add) {
        int absent;
	k = kh_put(s2i, d, ss, &absent);
        if (absent < 0) {
            return -1; // Out of memory
        } else if (absent) {
            char *ss_dup = strdup(ss);
            if (ss_dup) {
                kh_key(d, k) = ss_dup;
                kh_val(d, k) = kh_size(d) - 1;
            } else {
                kh_del(s2i, d, k);
		return -1; // Out of memory
            }
        }
    } else k = kh_get(s2i, d, ss);
    return k == kh_end(d)? -1 : kh_val(d, k);
}

int tbx_parse1(const tbx_conf_t *conf, int len, char *line, tbx_intv_t *intv)
{
    int i, b = 0, id = 1, ncols = 0;
    char *s;
    intv->ss = intv->se = 0; intv->beg = intv->end = -1;
    for (i = 0; i <= len; ++i) {
        if (line[i] == '\t' || line[i] == 0) {
            ++ncols;
            if (id == conf->sc) {
                intv->ss = line + b; intv->se = line + i;
            } else if (id == conf->bc) {
                // here ->beg is 0-based.
                intv->beg = intv->end = strtoll(line + b, &s, 0);
                if ( s==line+b ) return -1; // expected int
                if (!(conf->preset&TBX_UCSC)) --intv->beg;
                else ++intv->end;
                if (intv->beg < 0) intv->beg = 0;
                if (intv->end < 1) intv->end = 1;
            } else {
                if ((conf->preset&0xffff) == TBX_GENERIC) {
                    if (id == conf->ec)
                    {
                        intv->end = strtoll(line + b, &s, 0);
                        if ( s==line+b ) return -1; // expected int
                    }
                } else if ((conf->preset&0xffff) == TBX_SAM) {
                    if (id == 6) { // CIGAR
                        int l = 0;
                        char *t;
                        for (s = line + b; s < line + i;) {
                            long x = strtol(s, &t, 10);
                            char op = toupper_c(*t);
                            if (op == 'M' || op == 'D' || op == 'N') l += x;
                            s = t + 1;
                        }
                        if (l == 0) l = 1;
                        intv->end = intv->beg + l;
                    }
                } else if ((conf->preset&0xffff) == TBX_VCF) {
                    if (id == 4) {
                        if (b < i) intv->end = intv->beg + (i - b);
                    } else if (id == 8) { // look for "END="
                        int c = line[i];
                        line[i] = 0;
                        s = strstr(line + b, "END=");
                        if (s == line + b) s += 4;
                        else if (s) {
                            s = strstr(line + b, ";END=");
                            if (s) s += 5;
                        }
                        if (s && *s != '.') {
                            long long end = strtoll(s, &s, 0);
                            if (end <= intv->beg) {
                                static int reported = 0;
                                if (!reported) {
                                    int l = intv->ss ? (int) (intv->se - intv->ss) : 0;
                                    hts_log_warning("VCF INFO/END=%lld is smaller than POS at %.*s:%"PRIhts_pos"\n"
                                                    "This tag will be ignored. "
                                                    "Note: only one invalid END tag will be reported.",
                                                    end, l >= 0 ? l : 0,
                                                    intv->ss ? intv->ss : "",
                                                    intv->beg);
                                    reported = 1;
                                }
                            } else {
                                intv->end = end;
                            }
                        }
                        line[i] = c;
                    }
                }
            }
            b = i + 1;
            ++id;
        }
    }
    if (intv->ss == 0 || intv->se == 0 || intv->beg < 0 || intv->end < 0) return -1;
    return 0;
}

static inline int get_intv(tbx_t *tbx, kstring_t *str, tbx_intv_t *intv, int is_add)
{
    if (tbx_parse1(&tbx->conf, str->l, str->s, intv) == 0) {
        int c = *intv->se;
        *intv->se = '\0'; intv->tid = get_tid(tbx, intv->ss, is_add); *intv->se = c;
        if (intv->tid < 0) return -2;  // get_tid out of memory
        return (intv->beg >= 0 && intv->end >= 0)? 0 : -1;
    } else {
        char *type = NULL;
	switch (tbx->conf.preset&0xffff)
        {
            case TBX_SAM: type = "TBX_SAM"; break;
            case TBX_VCF: type = "TBX_VCF"; break;
            case TBX_UCSC: type = "TBX_UCSC"; break;
            default: type = "TBX_GENERIC"; break;
        }
        hts_log_error("Failed to parse %s, was wrong -p [type] used?\nThe offending line was: \"%s\"",
            type, str->s);
        return -1;
    }
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
    bidx_t *bidx = idx->bidx[i];
    lidx_t *lidx = &idx->lidx[i];
    khint_t k;
    int l;
    // the last entry is always valid
    for (l=lidx->n-2; l >= 0; l--) {
        if (lidx->offset[l] == (uint64_t)-1)
            lidx->offset[l] = lidx->offset[l+1];
    }
    if (bidx == 0) return;
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) // set loff
        if (kh_exist(bidx, k))
        {
            if ( kh_key(bidx, k) < idx->n_bins )
            {
                int bot_bin = hts_bin_bot(kh_key(bidx, k), idx->n_lvls);
                // disable linear index if bot_bin out of bounds
                kh_val(bidx, k).loff = bot_bin < lidx->n ? lidx->offset[bot_bin] : 0;
            }
            else
                kh_val(bidx, k).loff = 0;
        }
    if (free_lidx) {
        free(lidx->offset);
        lidx->m = lidx->n = 0;
        lidx->offset = 0;
    }
}

#define META_BIN(idx) ((idx)->n_bins + 1)

static inline int insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
    khint_t k;
    bins_t *l;
    int absent;
    k = kh_put(bin, b, bin, &absent);
    if (absent < 0) return -1; // Out of memory
    l = &kh_value(b, k);
    if (absent) {
	l->m = 1; l->n = 0;
        l->list = (hts_pair64_t*)calloc(l->m, sizeof(hts_pair64_t));
	if (!l->list) {
	    kh_del(bin, b, k);
            return -1;
        }
    } else if (l->n == l->m) {
        uint32_t new_m = l->m ? l->m << 1 : 1;
        hts_pair64_t *new_list = realloc(l->list, new_m * sizeof(hts_pair64_t));
        if (!new_list) return -1;
        l->list = new_list;
        l->m = new_m;
    }
    l->list[l->n].u = beg;
    l->list[l->n++].v = end;
    return 0;
}

int hts_idx_finish_chunk(hts_idx_t *idx, uint64_t final_offset)
{
    int i, ret = 0;
    if (idx == NULL || idx->z.finished) return 0; // do not run this function on an empty index or multiple times
    if (idx->z.save_tid >= 0) {
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, final_offset);
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
    }
    for (i = 0; i < idx->n; ++i) {
        update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
    }
    idx->z.finished = 1;

    return ret;
}

#endif /* TBX_EXTRA_H */
