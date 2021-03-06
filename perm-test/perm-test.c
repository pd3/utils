/* 
    Copyright (C) 2018-2020 Genome Research Ltd.
    
    Author: Petr Danecek <pd3@sanger.ac.uk>
    
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
    DEALINGS IN THE SOFTWARE.
*/
/*
    Permutation test. On input we are given a set of calls (regions of
    arbitrary length), a set of background regions (no enrichment expected),
    set of target regions (regions which are tested for enrichment), and
    the chromosome lengths.

    We randomly place the calls onto the regions (taking their length into
    account, therefore not just permuting labels) and calculate how many times
    the number of on-target-calls matches or exceeds the number of on-target-calls
    observed in input data.

    Because the regions can be small relative to the whole genome (imagine exomes
    which cover just 2% of the genome), we do not want to be placing the calls 
    genome wide, that would be a huge waste of CPU time. Instead, we create
    artificial chromosomes which condense the uninformative regions. Any placement
    of a call on this artificial genome hits a background or target region.
    The program correctly handles boundaries, both region and chromosomal.

    Because the artificial genome depends on the call length, it would become
    prohibitive to keep the genomes for all call lengths in memory at the same
    time with many calls. Therefore we iterate each call separately, one at a
    time, and instead accumulate the number of hits for each iteration.
*/

#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <regidx.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/khash_str2int.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <getopt.h>
#include "../libs/version.h"
#include "../libs/dist.h"

KHASH_MAP_INIT_INT(int2int, int)

#define NBIN2CHR (1<<16)

typedef struct
{
    uint32_t beg, len:31, is_tgt:1;
}
reg_t;

typedef struct
{
    char *name;     // chrom name
    uint32_t len;   // chrom length
    uint32_t clen;  // call length used to build the artificial chrom
    uint32_t alen;  // length of the artificial chrom
    uint32_t amax;  // the maximum coordinate on achrom before the call goes beyond the end
    reg_t *regs;    // list of spliced and merged regions
    uint32_t nregs;
    regidx_t *tgt_idx;  // index to the artificial chrom
    regidx_t *bg_idx;   // for the --hit-no-bg option
}
chr_t;

typedef struct
{
    uint32_t nobs_tgt_hits;     // number of observed hits in target regions
    uint32_t *niter_hits;       // number of hits in each permutation
    uint32_t *calls;            // call lengths, sorted in ascending order
    uint32_t max_call_len;      // skip calls longer than this
    uint32_t nbin2chr;
    chr_t *chr, **bin2chr;
    int nchr, ncalls;
    char *calls_fname, *background_fname, *targets_fname, *ref_fai_fname, *output_fname;
    uint32_t niter, nrounds;    // number of iterations to run and the number of rounds
    regidx_t *tgt_idx, *bg_idx;
    regitr_t *tgt_itr, *bg_itr;
    int debug, hit_no_bg, print_placements, ndist;
    FILE *out_fh;
}
args_t;

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

#define PUSH_REGION(_beg,_end,_is_tgt) \
{\
    if ( rep_end1 && (_beg) < rep_end1 ) (_beg) = rep_end1; \
    if ( (_beg) <= (_end) ) \
    { \
        nregs++; \
        hts_expand(reg_t,nregs,mregs,regs); \
        regs[nregs-1].beg = (_beg); \
        regs[nregs-1].len = (_end) - (_beg) + 1; \
        regs[nregs-1].is_tgt = (_is_tgt); \
        if ( rep_end1 < (_end) + 1 ) rep_end1 = (_end) + 1; \
        if ( (_is_tgt) && rep_tgt_end1 < (_end) + 1 ) rep_tgt_end1 = (_end) + 1; \
    } \
}

void merge_and_splice_regions(args_t *args, chr_t *chr)
{
    reg_t *regs = NULL;
    int i, nregs = 0, mregs = 0;

    uint32_t rep_tgt_end1 = 0; // the end of the last reported tgt region, 1-based real-chrom coordinates
    uint32_t rep_end1  = 0;    // the end of the last reported bg or tgt region, 1-based real-chrom coordinates; rep_end1 >= rep_tgt_end1

    // iterate over all bg regions, condense and splice as necessary
    if ( regidx_overlap(args->bg_idx,chr->name,0,REGIDX_MAX,args->bg_itr) )
    {
        while ( regitr_overlap(args->bg_itr) )
        {
            uint32_t bg_beg = args->bg_itr->beg;    // current bg region's 0-based real-chrom coordinates
            uint32_t bg_end = args->bg_itr->end;

            // adjust the start if the region overlaps a previously reported one
            if ( rep_end1 && rep_end1 - 1 >= bg_beg ) bg_beg = rep_end1;

            // check if there is a tgt region in the gap between this bg and the previously reported region
            if ( rep_end1 < bg_beg && regidx_overlap(args->tgt_idx,chr->name,rep_end1,bg_beg-1,args->tgt_itr) )
            {
                while ( regitr_overlap(args->tgt_itr) )
                {
                    PUSH_REGION(args->tgt_itr->beg, args->tgt_itr->end, 1);
                }

                if ( rep_end1 && rep_end1 - 1 >= bg_beg ) bg_beg = rep_end1;
            }
            if ( bg_beg > bg_end ) continue;    // the bg region was overwritten by tgt regions completely

            // check if there is a tgt region which overlaps this bg region (or what's left of it)
            if ( regidx_overlap(args->tgt_idx,chr->name,bg_beg,bg_end,args->tgt_itr) )
            {
                while ( bg_beg <= bg_end && regitr_overlap(args->tgt_itr) )
                {
                    // add the first part of the bg region
                    PUSH_REGION(bg_beg, args->tgt_itr->beg - 1, 0);

                    // and the tgt region itself
                    PUSH_REGION(args->tgt_itr->beg, args->tgt_itr->end, 1);

                    bg_beg = rep_end1;
                }
            }

            // flush all what is left of the bg region
            PUSH_REGION(bg_beg, bg_end, 0);
        }
    }

    // flush all target regions which start beyond the background regions
    if ( regidx_overlap(args->tgt_idx,chr->name,rep_end1,REGIDX_MAX, args->tgt_itr) )
    {
        while ( regitr_overlap(args->tgt_itr) )
            PUSH_REGION(args->tgt_itr->beg, args->tgt_itr->end, 1);
    }

    // merge consecutive regions into one
    for (i=1; i<nregs; i++)
    {
        if ( regs[i-1].is_tgt != regs[i].is_tgt ) continue;
        if ( regs[i-1].beg + regs[i-1].len < regs[i].beg ) continue;
        assert( regs[i-1].beg + regs[i-1].len - 1 < regs[i].beg + regs[i].len - 1 );
        regs[i-1].len = regs[i].beg + regs[i].len - regs[i-1].beg;
        if ( i + 1 < nregs )
            memmove(&regs[i], &regs[i+1], (nregs-i-1)*sizeof(*regs));
        nregs--; i--;
    }

    // downsize
    chr->nregs = nregs;
    chr->regs  = (reg_t*) realloc(regs, nregs*sizeof(*chr->regs));

    if ( args->debug )
    {
        for (i=0; i<chr->nregs; i++)
            fprintf(args->out_fh, "%s\t%s\t%d\t%d\n", chr->regs[i].is_tgt ? "TGT" : "BG", chr->name, chr->regs[i].beg+1, chr->regs[i].beg + chr->regs[i].len);
    }
}
static void init_bin2chr(args_t *args)
{
    args->bin2chr = (chr_t**) calloc(NBIN2CHR,sizeof(*args->bin2chr));
    double tot_len = 0;
    int i,j;
    for (i=0; i<args->nchr; i++)
        tot_len += args->chr[i].len;

    double bin_len = tot_len / NBIN2CHR, max_diff = 0;
    for (i=0,j=0; i<args->nchr; i++)
    {
        int k, kmax = args->chr[i].len/tot_len*NBIN2CHR;
        for (k=0; k<kmax; k++)
            args->bin2chr[j++] = &args->chr[i];
        double diff = fabs(kmax*bin_len - args->chr[i].len)*100./args->chr[i].len;
        if ( max_diff < diff ) max_diff = diff;
    }
    args->nbin2chr = j;
    fprintf(args->out_fh,"MSG\tMaximum chromosome randomization error due to length discretization: %.1e%%\n", max_diff);
}
int is_bed_file(const char *fname)
{
    int len = strlen(fname);
    if ( len>=7 && !strcasecmp(".bed.gz",fname+len-7) ) return 1; 
    else if ( len>=8 && !strcasecmp(".bed.bgz",fname+len-8) ) return 1;
    else if ( len>=4 && !strcasecmp(".bed",fname+len-4) ) return 1;
    return 0;
}

// this is to sort the calls by size
static int uint32t_cmp(const void *aptr, const void *bptr)
{
    uint32_t a = *((uint32_t*)aptr);
    uint32_t b = *((uint32_t*)bptr);
    if ( a < b ) return -1;
    if ( a > b ) return 1;
    return 0;
}
void init_data(args_t *args)
{
    args->niter_hits = (uint32_t*) calloc(args->niter, sizeof(*args->niter_hits));
    if ( !args->niter_hits )
        error("Could not allocate %zu bytes. Consider using the second argument to --niter.\n", args->niter*sizeof(*args->niter_hits));

    args->tgt_idx = regidx_init(args->targets_fname,NULL,NULL,0,NULL);
    if ( !args->tgt_idx )
        error("Error indexing %s\n", args->targets_fname);
    args->bg_idx  = regidx_init(args->background_fname,NULL,NULL,0,NULL);
    if ( !args->bg_idx )
        error("Error indexing %s\n", args->background_fname);
    args->tgt_itr = regitr_init(args->tgt_idx);
    args->bg_itr  = regitr_init(args->bg_idx);

    // get a list of all chromosomes from both tgt and bg regions
    void *chr_hash = khash_str2int_init();
    int i,n;
    char **names = regidx_seq_names(args->bg_idx, &args->nchr);
    args->chr = (chr_t*) calloc(args->nchr,sizeof(*args->chr));
    for (i=0; i<args->nchr; i++)
    {
        args->chr[i].name = strdup(names[i]);
        khash_str2int_set(chr_hash, args->chr[i].name, i);
    }
    names = regidx_seq_names(args->tgt_idx, &n);
    for (i=0; i<n; i++)
    {
        if ( khash_str2int_has_key(chr_hash,names[i]) ) continue;
        args->nchr++;
        args->chr = (chr_t*) realloc(args->chr,sizeof(*args->chr)*args->nchr);
        chr_t *chr = &args->chr[args->nchr-1];
        memset(chr, 0, sizeof(chr_t));
        chr->name = strdup(names[i]);
        khash_str2int_set(chr_hash, chr->name, args->nchr-1);
    }

    // merge and splice bg and tgt regions
    for (i=0; i<args->nchr; i++)
        merge_and_splice_regions(args, &args->chr[i]);

    // read chromosome lengths
    BGZF *fp = bgzf_open(args->ref_fai_fname,"r");
    if ( !fp ) error("Failed to read: %s\n", args->ref_fai_fname);
    kstring_t str = {0,0,0};
    while ( bgzf_getline(fp, '\n', &str) > 0 )
    {
        char keep, *tmp, *ptr = str.s;
        while ( *ptr && !isspace(*ptr) ) ptr++;
        if ( !*ptr ) error("Could not parse %s: %s\n", args->ref_fai_fname,str.s);
        keep = *ptr; *ptr = 0;
        int ichr;
        if ( khash_str2int_get(chr_hash,str.s,&ichr)!=0 ) continue;
        *ptr = keep;
        chr_t *chr = &args->chr[ichr];
        ptr++;
        while ( *ptr && isspace(*ptr) ) ptr++;
        chr->len = strtol(ptr, &tmp, 10);
        if ( tmp==ptr ) error("Could not parse %s: %s\n", args->ref_fai_fname,str.s);
    }
    if ( bgzf_close(fp)!=0 ) error("close failed: %s\n", args->ref_fai_fname);
    khash_str2int_destroy(chr_hash);

    for (i=0; i<args->nchr; i++)
    {
        if ( args->chr[i].len ) continue;
        fprintf(stderr,"Warning: Skipping the sequence \"%s\" because it is not listed in %s\n",args->chr[i].name,args->ref_fai_fname);
        chr_t *chr = &args->chr[i];
        free(chr->name);
        free(chr->regs);
        if ( args->nchr - i > 1 ) memmove(&args->chr[i],&args->chr[i+1],sizeof(*args->chr)*(args->nchr - i - 1));
        args->nchr--;
    }
    init_bin2chr(args);

    // read calls
    int nskip  = 0;
    int is_bed = is_bed_file(args->calls_fname);
    fp = bgzf_open(args->calls_fname,"r");
    if ( !fp ) error("Failed to read: %s\n", args->calls_fname);
    while ( bgzf_getline(fp, '\n', &str) > 0 )
    {
        char *chr_beg, *chr_end;
        uint32_t beg, end;
        int ret = regidx_parse_tab(str.s, &chr_beg, &chr_end, &beg, &end, NULL, NULL);
        if ( ret < -1 ) error("Could not parse %s: %s\n", args->calls_fname, str.s);
        if ( is_bed ) beg++;
        chr_end[1] = 0;

        // Check if the call overlaps a bg or tgt region
        int is_tgt = regidx_overlap(args->tgt_idx, chr_beg,beg,end, NULL);
        int is_bg  = regidx_overlap(args->bg_idx, chr_beg,beg,end, NULL);
        if ( (!is_tgt && !is_bg) || end - beg + 1 > args->max_call_len )
        {
            nskip++;
            if ( args->debug ) fprintf(args->out_fh, "CALL\tSKIP\t%s\t%d\t%d\n", chr_beg,beg+1,end+1);
            continue;
        }
        if ( is_tgt && (!args->hit_no_bg || !is_bg) ) args->nobs_tgt_hits++;
        if ( args->debug )
            fprintf(args->out_fh, "CALL\t%s\t%s\t%d\t%d\n",is_tgt && is_bg ? "TGT_BG" : ( is_tgt ? "TGT" : "BG"),chr_beg,beg+1,end+1);

        args->ncalls++;
        args->calls = (uint32_t*) realloc(args->calls, sizeof(*args->calls)*args->ncalls);
        if ( !args->calls ) error("Could not alloc %zu bytes\n", sizeof(*args->calls)*args->ncalls);
        args->calls[args->ncalls-1] = end - beg + 1;
    }
    if ( bgzf_close(fp)!=0 ) error("close failed: %s\n", args->calls_fname);
    free(str.s);

    if ( !args->ncalls && !args->debug )
        error("Error: none of the calls intersects the tgt or bg regions (calls larger than `-m %d` were excluded)\n", args->max_call_len);
    qsort(args->calls, args->ncalls, sizeof(*args->calls), uint32t_cmp);

    fprintf(args->out_fh, "NCALLS\t%d\t%d\n", args->ncalls,nskip);
}
void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nchr; i++)
    {
        chr_t *chr = &args->chr[i];
        free(chr->name);
        free(chr->regs);
        if ( chr->tgt_idx ) regidx_destroy(chr->tgt_idx);
        if ( chr->bg_idx ) regidx_destroy(chr->bg_idx);
    }
    free(args->chr);
    free(args->niter_hits);
    free(args->calls);
    free(args->bin2chr);
    regidx_destroy(args->tgt_idx);
    regidx_destroy(args->bg_idx);
    regitr_destroy(args->tgt_itr);
    regitr_destroy(args->bg_itr);
    if ( args->out_fh && fclose(args->out_fh)!=0 )
        error("Error: close failed on %s\n", args->output_fname ? args->output_fname : "standard outpout");
}
static inline void update_alen_amax(chr_t *chr, uint32_t reg_len, uint32_t reg_beg0)
{
    // Check if there is enough space for the call when placed within the
    // half-open interval [reg_beg0,reg_beg0+reg_len) until the end of the chromosome
    // and set rand_max if not (coordinates are 0-based)
    if ( !chr->amax && chr->len < reg_beg0 + reg_len + chr->clen - 1 )
    {
        if ( chr->len >= reg_beg0 + chr->clen )
            chr->amax = chr->alen + chr->len - chr->clen - reg_beg0;    // amax is 0-based
        else
            chr->amax = chr->alen - 1;
    }
    chr->alen += reg_len;
}
void init_chr(args_t *args, chr_t *chr, uint32_t call_len)
{
    if ( call_len <=0 ) error("Error: call length %d??\n", call_len);

    chr->amax = 0;
    chr->alen = 0;
    chr->clen = call_len;

    size_t payload_size = args->print_placements ? sizeof(uint32_t) : 0;

    if ( chr->tgt_idx ) regidx_destroy(chr->tgt_idx);
    chr->tgt_idx = regidx_init(NULL,NULL,NULL,payload_size,NULL);
    if ( args->hit_no_bg || args->print_placements )
    {
        if ( chr->bg_idx ) regidx_destroy(chr->bg_idx);
        chr->bg_idx = regidx_init(NULL,NULL,NULL,payload_size,NULL);
    }

    uint32_t i, rep_end1 = 0, clen1 = call_len - 1;
    char *chr_name_end = chr->name + strlen(chr->name) - 1;

    if ( clen1 == 0 )
    {
        for (i=0; i<chr->nregs; i++)
        {
            if ( args->debug > 1 )
                fprintf(args->out_fh, "L%d_%s\t%s\t%d\t%d\t..\t%d\t%d\n", call_len, chr->regs[i].is_tgt ? "TGT" : "BG", chr->name, chr->regs[i].beg+1, chr->regs[i].beg + chr->regs[i].len, chr->alen + 1, chr->alen + chr->regs[i].len);

            void *payload = args->print_placements ? &chr->regs[i].beg : NULL;

            if ( chr->regs[i].is_tgt )
                regidx_push(chr->tgt_idx, chr->name, chr_name_end, chr->alen, chr->alen + chr->regs[i].len - 1, payload);
            else if ( args->hit_no_bg || args->print_placements )
                regidx_push(chr->bg_idx, chr->name, chr_name_end, chr->alen, chr->alen + chr->regs[i].len - 1, payload);

            chr->alen += chr->regs[i].len;
        }
        chr->amax = chr->alen - 1;
    }
    else
    {
        for (i=0; i<chr->nregs; i++)
        {
            // the left overhang can be clen1 bp or shorter if the call is bigger than the gap between two regions
            uint32_t pbeg = chr->regs[i].beg >= clen1 ? chr->regs[i].beg - clen1 : 0;
            if ( pbeg < rep_end1 ) pbeg = rep_end1;
            if ( pbeg < chr->regs[i].beg )
                update_alen_amax(chr, chr->regs[i].beg - pbeg, pbeg);

            if ( args->debug > 2 )
                fprintf(args->out_fh,"L%d_%s\t%s\t%d\t%d\t..\t%d\t%d\n", call_len, chr->regs[i].is_tgt ? "TGT" : "BG", chr->name, chr->regs[i].beg+1, chr->regs[i].beg + chr->regs[i].len, chr->alen + 1, chr->alen + chr->regs[i].len);

            void *payload = args->print_placements ? &chr->regs[i].beg : NULL;

            if ( chr->regs[i].is_tgt )
                regidx_push(chr->tgt_idx, chr->name, chr_name_end, chr->alen, chr->alen + chr->regs[i].len - 1, payload);
            else if ( args->hit_no_bg || args->print_placements )
                regidx_push(chr->bg_idx, chr->name, chr_name_end, chr->alen, chr->alen + chr->regs[i].len - 1, payload);

            update_alen_amax(chr, chr->regs[i].len, chr->regs[i].beg);
            rep_end1 = chr->regs[i].beg + chr->regs[i].len;
        }
        if ( !chr->amax ) chr->amax = chr->alen - 1;
    }
    if ( args->debug > 2 )
        fprintf(args->out_fh, "ALEN_AMAX%d\t%d\t%d\n", call_len, chr->alen,chr->amax);
}
void run_test(args_t *args, uint32_t call_len)
{
    regitr_t *tgt_itr = args->print_placements ? args->tgt_itr : NULL;
    regitr_t *bg_itr  = args->print_placements ? args->bg_itr  : NULL;

    int i;
    for (i=0; i<args->niter; i++)
    {
        // randomly assign a chromosome
        int ibin = (double)random()/((double)RAND_MAX+1) * args->nbin2chr;
        chr_t *chr = args->bin2chr[ibin];

        if ( chr->len <= call_len )
        {
            if ( args->debug > 3 ) fprintf(args->out_fh, "RAND\t%s\t%u\t%u\n", chr->name,1,chr->len);
            if ( regidx_overlap(args->tgt_idx,chr->name,0,chr->len,tgt_itr) )
            {
                if ( args->print_placements && regitr_overlap(tgt_itr) )
                    fprintf(args->out_fh, "POS\t%s\t%d\t%d\t1\n", chr->name,1,call_len);
                if ( !args->hit_no_bg || !regidx_overlap(chr->bg_idx,chr->name,0,chr->len,NULL) ) args->niter_hits[i]++;
            }
            else if ( args->print_placements && regidx_overlap(args->bg_idx,chr->name,0,chr->len,bg_itr) && regitr_overlap(bg_itr) )
                fprintf(args->out_fh, "POS\t%s\t%d\t%d\t0\n", chr->name,1,call_len);
            continue;
        }

        if ( !chr->tgt_idx || chr->clen!=call_len ) init_chr(args, chr, call_len);

        uint32_t pos = (double)random()/((double)RAND_MAX+1) * (chr->amax + 1);
        if ( args->debug > 3 ) fprintf(args->out_fh, "RAND\t%s\t%u\t%u\n", chr->name,chr->amax+1,pos+1);
        if ( regidx_overlap(chr->tgt_idx, chr->name, pos, pos + call_len - 1, tgt_itr) )
        {
            if ( args->print_placements && regitr_overlap(tgt_itr) )
            {
                uint32_t ori_pos = regitr_payload(tgt_itr,uint32_t) - tgt_itr->beg + pos + 1;
                fprintf(args->out_fh, "POS\t%s\t%d\t%d\t1\n", chr->name,ori_pos,ori_pos+call_len-1);
            }
            if ( !args->hit_no_bg || !regidx_overlap(chr->bg_idx,chr->name,pos,pos+call_len-1,NULL) ) args->niter_hits[i]++;
        }
        else if ( args->print_placements && regidx_overlap(chr->bg_idx,chr->name,pos,pos+call_len-1,bg_itr) && regitr_overlap(bg_itr) )
        {
            uint32_t ori_pos = regitr_payload(bg_itr,uint32_t) - bg_itr->beg + pos + 1;
            fprintf(args->out_fh, "POS\t%s\t%d\t%d\t0\n", chr->name,ori_pos,ori_pos+call_len-1);
        }
    }
}

static void usage(void)
{
    error(
        "\n"
        "Program: perm-test %s\n"
        "License: The MIT/Expat license\n"
        "This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        "About: Run enrichment (permutation) test. It is efficient, does not attempt to place calls\n"
        "       in inaccessible regions, and correctly handles calls and regions of different sizes.\n"
        "       All coordinates are 1-based, inclusive, unless the file names have \".bed\", \".bed.gz\"\n"
        "       or \".bed.bgz\" extension.\n"
        "Usage: perm-test [OPTIONS]\n"
        "Options:\n"
        "   -b, --background-regs FILE      Background regions, expected not to be enriched: chr,beg,end\n"
        "   -c, --calls FILE                Calls: chr,beg,end\n"
        "   -d, --debug-regions             Print the spliced regions (and stop)\n"
        "   -f, --ref-fai FILE              Chromosome lengths, given for example as fai index: chr,length\n"
        "   -m, --max-call-length INT       Skip big calls [10e6]\n"
        "   -n, --niter NUM1[,NUM2]         Number of iterations: total (NUM1) and per-batch (NUM2, smaller number reduces memory)\n"
        "       --no-bg-overlap             Permuted variants must not overlap background regions\n"
        "       --nprecise-dist INT         Orders of magnitude to represent exactly in the distribution [5]\n"
        "   -o, --output FILE               Place output in FILE\n"
        "   -p, --print-placements          Print all random placements\n"
        "   -s, --random-seed INT           Random seed\n"
        "   -t, --target-regs FILE          Target regions, expected to be enriched: chr,beg,end\n"
        "\n"
        "Examples:\n"
        "   # Take a set of \"background\" regions (e.g. exome baits), \"target\" regions (e.g. selected genes)\n"
        "   # and calculate how likely it is for a set of calls (e.g. CNVs) to overlap the target regions given\n"
        "   # the accessible background+target regions. The program ensures that overlapping and duplicate\n"
        "   # regions are handled correctly.\n"
        "   perm-test -b baits.txt -t genes.txt -c calls.txt -f ref.fai -n 1e9,1e8\n"
        "\n"
        "   # Same as above but this time count only permuted variants that overlap target regions but do not\n"
        "   # overlap any background region (e.g. CNVs that hit only a non-coding region but not a coding sequence).\n"
        "   # Note that target regions take precedence over background regions in case the sets overlap. If that is\n"
        "   # not desired, use `bedtools subtract -a A.bed -b B.bed` to clip and write regions from file A not present\n"
        "   # in B (e.g. `bedtools substract -a baits.bed -b coding.bed > noncoding.bed`).\n"
        "   perm-test -b coding.bed -t noncoding.bed -c calls.txt -f ref.fai -n 1e9,1e8 --no-bg-overlap\n"
        "\n",
        UTILS_VERSION
        );
}

int main(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->nrounds = 1;
    args->max_call_len = 10e6;
    args->ndist = 5;
    static struct option loptions[] =
    {
        {"nprecise-dist",required_argument,NULL,2},
        {"print-placements",no_argument,NULL,'p'},
        {"no-bg-overlap",no_argument,NULL,1},
        {"max-call-len",required_argument,NULL,'m'},
        {"ref-fai",required_argument,NULL,'f'},
        {"output",required_argument,NULL,'o'},
        {"ninter",required_argument,NULL,'n'},
        {"debug-regions",required_argument,NULL,'d'},
        {"calls",required_argument,NULL,'c'},
        {"background-regs",required_argument,NULL,'b'},
        {"target-regs",required_argument,NULL,'t'},
        {"random-seed",required_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int c, seed = -1;
    while ((c = getopt_long(argc, argv, "?hc:b:t:n:ds:f:o:m:p",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 : args->hit_no_bg = 1; break;
            case  2 : 
                args->ndist = strtol(optarg, &tmp, 10);
                if ( tmp==optarg || *tmp ) error("Could not parse --nprecise-dist %s\n", optarg);
                if ( args->ndist<1 ) error("Nonsense value to --nprecise-dist %s\n", optarg);
                break;
            case 'p': args->print_placements = 1; break;
            case 'm': 
                args->max_call_len = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --max-call-length %s\n", optarg);
                break;
            case 's': 
                seed = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --random-seed %s\n", optarg);
                break;
            case 'n': 
                {
                    uint64_t nbatch = 0, ntot_iter = strtod(optarg, &tmp); 
                    if ( tmp!=optarg && *tmp==',' )
                    {
                        char *tmp2 = tmp+1;
                        nbatch = strtod(tmp2, &tmp);
                        if ( tmp==tmp2 || *tmp ) error("Could not parse --niter %s\n", optarg);
                    }
                    else if ( tmp==optarg ) error("Could not parse --niter %s\n", optarg);
                    args->niter = nbatch ? nbatch : ntot_iter;
                    args->nrounds = ntot_iter % args->niter ? ntot_iter / args->niter + 1 : ntot_iter / args->niter;
                    break;
                }
            case 'f': args->ref_fai_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'c': args->calls_fname = optarg; break;
            case 't': args->targets_fname = optarg; break;
            case 'b': args->background_fname = optarg; break;
            case 'd': args->debug++; break;
            case 'h':
            case '?':
            default: usage(); break;
        }
    }
    if ( argc < 2 ) usage();
    if ( !args->ref_fai_fname ) error("Missing the -f option\n");
    if ( !args->calls_fname ) error("Missing the -c option\n");
    if ( !args->targets_fname ) error("Missing the -t option\n");
    if ( !args->background_fname ) error("Missing the -b option\n");
    args->out_fh = args->output_fname ? fopen(args->output_fname,"w") : stdout;
    if ( !args->out_fh ) error("Error: Could not open %s for writing\n", args->output_fname);
    if ( seed==-1 )
    {
        struct timeval time; 
        gettimeofday(&time,NULL);
        seed = time.tv_sec + time.tv_usec;
    }
    if ( args->niter )
    {
        fprintf(stderr,"Using random seed %d\nRunning %u rounds, %e iterations each\n", seed, args->nrounds, (double)args->niter);
        fprintf(args->out_fh, "# NITER_ROUNDS:\n");
        fprintf(args->out_fh, "#    - number of iterations total\n");
        fprintf(args->out_fh, "#    - number of batches\n");
        fprintf(args->out_fh, "# NCALLS:\n");
        fprintf(args->out_fh, "#    - number of calls used\n");
        fprintf(args->out_fh, "#    - number of calls skipped because of their size or inaccessibility\n");
        fprintf(args->out_fh, "# TEST_ENR:\n");
        fprintf(args->out_fh, "#    - number of iterations\n");
        fprintf(args->out_fh, "#    - number of times the simulation had the same or more hits than observed in input data\n");
        fprintf(args->out_fh, "#    - P-value derived from the two fields\n");
        fprintf(args->out_fh, "# TEST_DPL:\n");
        fprintf(args->out_fh, "#    - number of iterations\n");
        fprintf(args->out_fh, "#    - number of times the simulation had the same or fewer hits than observed in input data\n");
        fprintf(args->out_fh, "#    - P-value derived from the two fields\n");
        fprintf(args->out_fh, "# TEST_FOLD:\n");
        fprintf(args->out_fh, "#    - number of hits in the input data\n");
        fprintf(args->out_fh, "#    - average number of hits in simulations\n");
        fprintf(args->out_fh, "#    - average of standard deviations approximated from each batch using the current average estimate\n");
        fprintf(args->out_fh, "# DIST, the distribution of number of hits in the simulations:\n");
        fprintf(args->out_fh, "#    - numbef of hits in a simulation: the start of the half-closed, half-open distribution bin [beg,end)\n");
        fprintf(args->out_fh, "#    - numbef of hits in a simulation: the end of the half-closed, half-open distributin bin [beg,end)\n");
        fprintf(args->out_fh, "#    - number of simulations with this many hits (raw count)\n");
        fprintf(args->out_fh, "#    - number of simulations with this many hits (density)\n");
        fprintf(args->out_fh, "VERSION\t%s\n",UTILS_VERSION);
        fprintf(args->out_fh, "CMD\t%s",argv[0]); for (c=1; c<argc; c++) fprintf(args->out_fh, " %s", argv[c]); fprintf(args->out_fh, "\n");
        fprintf(args->out_fh, "SEED\t%d\n", seed);
        fprintf(args->out_fh, "NITER_ROUNDS\t%e\t%u\n", (double)args->niter,args->nrounds);
    }
    srandom(seed);

    init_data(args);

    if ( args->debug!=1 )
    {
        dist_t *dist = dist_init(args->ndist);
        double navg_tgt_hits = 0, dev = 0;
        uint64_t nexc = 0, nfew = 0, ntot = 0;
        uint32_t i,j;
        for (j=0; j<args->nrounds; j++)
        {
            memset(args->niter_hits, 0, sizeof(*args->niter_hits)*args->niter);
            for (i=0; i<args->ncalls; i++)
                run_test(args, args->calls[i]);

            if ( !args->niter ) break;  // this is can be true only when debugging

            for (i=0; i<args->niter; i++)
            {
                if ( args->niter_hits[i] >= args->nobs_tgt_hits ) nexc++;
                if ( args->niter_hits[i] <= args->nobs_tgt_hits ) nfew++;
                navg_tgt_hits += args->niter_hits[i];
            }
            ntot += args->niter;

            double avg = navg_tgt_hits/ntot;
            for (i=0; i<args->niter; i++)
                dev += (args->niter_hits[i] - avg)*(args->niter_hits[i] - avg);

            for (i=0; i<args->niter; i++)
                dist_insert(dist, args->niter_hits[i]);
        }
        if ( ntot )
        {
            double pval_enr = nexc ? (double)nexc/ntot : (double)1/ntot;
            double pval_dpl = nfew ? (double)nfew/ntot : (double)1/ntot;
            fprintf(args->out_fh, "TEST_ENR\t%"PRIu64"\t%"PRIu64"\t%s%e\n", ntot,nexc,nexc?"":"<",pval_enr);
            fprintf(args->out_fh, "TEST_DPL\t%"PRIu64"\t%"PRIu64"\t%s%e\n", ntot,nfew,nfew?"":"<",pval_dpl);
            fprintf(args->out_fh, "TEST_FOLD\t%"PRIu32"\t%f\t%f\n", args->nobs_tgt_hits,navg_tgt_hits/ntot,sqrt(dev/ntot));

            int n = dist_n(dist);
            uint32_t beg, end;
            for (i=0; i<n; i++)
            {
                uint64_t cnt = dist_get(dist, i, &beg, &end);
                if ( !cnt ) continue;
                fprintf(args->out_fh, "DIST\t%u\t%u\t%"PRIu64"\t%e\n", beg, end, cnt, (double)cnt/(end-beg));
            }
        }
        dist_destroy(dist);
    }

    destroy_data(args);
    free(args);

    return 0;
}

