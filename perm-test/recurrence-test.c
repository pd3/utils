/* 
    Copyright (C) 2020 Genome Research Ltd.
    
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
    Recurrence test.
    
    On input we are given a set of calls (regions of arbitrary length), a set of
    labeled regions (e.g. exons in a set of genes), a set of accessible regions
    (e.g. exome baits), and chromosome lengths.

    In each iteration all calls are randomly placed on the genome and for each
    label (gene) we count how many times it was hit by more/fewer calls than
    observed. Note that we evaluate probability conditioned on all calls being
    placed in accessible regions: each call is placed repeatedly until it hits
    an accessible region, placements outside of accessible regions do not count.

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
#include <htslib/kbitset.h>
#include <getopt.h>
#include "../libs/version.h"

#define DIST_N(args) (((args)->ncall+1)*(args)->regs.nlbl)
#define DIST_IDX(args,ilbl,cnt) (((args)->ncall+1)*(ilbl)+(cnt))

typedef struct
{
    int ichr, ilbl;
    uint32_t beg,end;   // 0-based coordinates, inclusive
}
reg_t;

typedef struct
{
    char *name;     // chrom name
    uint32_t len;   // chrom length
}
chr_t;

typedef struct
{
    uint32_t nreg, mreg;
    uint32_t nlbl, nchr;
    reg_t *reg;
    chr_t *chr;
    char **lbl;
    void *lbl2i;
    void *chr2i;
}
regs_t;

typedef struct
{
    regs_t regs;

    uint32_t
        ncall, *call,           // stores only call lengths, nothing else is needed
        max_call_len,
        niter,
        ntry,                   // number of random placements into inaccessible regions before giving up
        *nobs, *neq,            // arrays of regs.nlbl gene counters: observed in data, equal to observed,
        *nexc, *nfew,           //  exceeded, or fewer than observed
        *nhit,                  // number of hits in one iteration
        *dist;                  // distribution of hits, see the DIST_* macros for indexing
    kbitset_t *hit;             // bitmask of genes hit by a single call
    void *gene2i;
    char **i2gene;

    regidx_t *tgt_idx, *bg_idx;
    regitr_t *tgt_itr;

    char *calls_fname, *accessible_fname, *labeled_fname, *ref_fai_fname, *output_fname;
    int debug;
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

static int is_bed_file(const char *fname)
{
    int len = strlen(fname);
    if ( len>=7 && !strcasecmp(".bed.gz",fname+len-7) ) return 1;
    else if ( len>=8 && !strcasecmp(".bed.bgz",fname+len-8) ) return 1;
    else if ( len>=4 && !strcasecmp(".bed",fname+len-4) ) return 1;
    return 0;
}
static int parse_line(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, char **lbl_beg, char **lbl_end)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments
    
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se ) { fprintf(stderr,"Could not parse the line: %s\n", line); return -2; }

    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) { fprintf(stderr,"Could not parse the line: %s\n", line); return -2; }
    if ( *beg==0 ) { fprintf(stderr,"Could not parse the line, expected 1-based coordinate: %s\n", line); return -2; }
    (*beg)--;

    if ( !se[0] || !se[1] ) { fprintf(stderr,"Could not parse the line: %s\n", line); return -2; }
    ss = se+1;
    *end = strtod(ss, &se);
    if ( ss==se || (*se && !isspace(*se)) ) *end = *beg;
    else if ( *end==0 ) { fprintf(stderr,"Could not parse the line, expected 1-based coordinate: %s\n", line); return -2; }
    else (*end)--;

    if ( lbl_beg )
    {
        if ( !*se ) { fprintf(stderr,"Could not parse the line, expected a label: %s\n", line); return -2; }
        while ( *se && isspace(*se) ) se++;
        *lbl_beg = se;
        while ( *se && !isspace(*se) ) se++;
        *lbl_end = *se ? se : se-1;
        if ( *lbl_beg > *lbl_end ) { fprintf(stderr,"Could not parse the line, expected a label: %s\n", line); return -2; }
    }
    return 0;
}
static void read_regions(char *fname, regs_t *regs)
{
    void *chr2i = regs->chr2i;
    void *lbl2i = regs->lbl2i;
    kstring_t str = {0,0,0};
    BGZF *fp = bgzf_open(fname,"r");
    if ( !fp ) error("Failed to read: %s\n", fname);
    int is_bed = is_bed_file(fname);
    while ( bgzf_getline(fp, '\n', &str) > 0 )
    {
        char *chr_beg, *chr_end;
        char *lbl_beg, *lbl_end;
        uint32_t beg, end;
        int ret = parse_line(str.s, &chr_beg, &chr_end, &beg, &end, lbl2i ? &lbl_beg : NULL, lbl2i ? &lbl_end : NULL);
        if ( ret < -1 ) error("Could not parse %s: %s\n", fname, str.s);
        if ( ret < 0 ) continue;
        if ( is_bed ) beg++;
        chr_end[1] = 0;

        int ichr;
        if ( khash_str2int_get(chr2i,str.s,&ichr)!=0 ) continue;  // not listed in .fai, ignore

        regs->nreg++;
        hts_expand(reg_t,regs->nreg,regs->mreg,regs->reg);
        reg_t *reg = &regs->reg[regs->nreg-1];
        reg->beg  = beg;
        reg->end  = end;
        reg->ichr = ichr;

        if ( lbl2i )
        {
            lbl_end[1] = 0;
            if ( khash_str2int_get(lbl2i,lbl_beg,&reg->ilbl)!=0 )
            {
                regs->nlbl++;
                regs->lbl = (char**)realloc(regs->lbl,sizeof(*regs->lbl)*regs->nlbl);
                regs->lbl[regs->nlbl-1] = strdup(lbl_beg);
                khash_str2int_inc(lbl2i, regs->lbl[regs->nlbl-1]);
            }
            if ( khash_str2int_get(lbl2i,lbl_beg,&reg->ilbl)!=0 ) error("Failed to hash \"%s\"??\n", lbl_beg);
        }
        else
            reg->ilbl = 0;
    }
    if ( bgzf_close(fp)!=0 ) error("close failed: %s\n", fname);
    free(str.s);
}
static int cmp_regs(const void *aptr, const void *bptr)
{
    reg_t *a = (reg_t*)aptr;
    reg_t *b = (reg_t*)bptr;
    if ( a->ilbl < b->ilbl ) return -1;
    if ( a->ilbl > b->ilbl ) return 1;
    if ( a->ichr < b->ichr ) return -1;
    if ( a->ichr > b->ichr ) return 1;
    if ( a->beg < b->beg ) return -1;
    if ( a->beg > b->beg ) return 1;
    if ( a->end < b->end ) return -1;
    if ( a->end > b->end ) return 1;
    return 0;
}
static void merge_regions(regs_t *regs)
{
    qsort(regs->reg, regs->nreg, sizeof(*regs->reg), cmp_regs);

    reg_t *prev = NULL;
    int i, iprev = -1;
    for (i=0; i<regs->nreg; i++)
    {
        reg_t *reg  = &regs->reg[i];
        if ( !prev )
        {
            iprev = i;
            prev  = reg;
            continue;
        }
        if ( prev->ilbl == reg->ilbl && prev->ichr == reg->ichr && prev->end + 1 >= reg->beg )
        {
            if ( prev->end < reg->end ) prev->end = reg->end;
            continue;
        }
        if ( i - iprev > 1 )
        {
            memmove(&regs->reg[iprev+1], &regs->reg[i], sizeof(*regs->reg)*(regs->nreg - i));
            regs->nreg -= i - iprev - 1;
            i = iprev;
        }
        prev = NULL;
    }
    if ( prev && iprev + 1 < regs->nreg ) regs->nreg = iprev + 1;
}
static void trim_chrs(regs_t *regs, regidx_t *idx)
{
    int i;
    for (i=0; i<regs->nchr; i++)
    {
        if ( regidx_overlap(idx,regs->chr[i].name,0,regs->chr[i].len-1,NULL) ) continue;
        free(regs->chr[i].name);
        if ( i+1 < regs->nchr ) memmove(&regs->chr[i],&regs->chr[i+1],sizeof(*regs->chr)*(regs->nchr-i-1));
        regs->nchr--;
    }
}
static void init_data(args_t *args)
{
    int i;
    uint32_t genome_len = 0, bg_len = 0;

    regs_t *regs = &args->regs;
    regs->chr2i = khash_str2int_init();
    
    // read .fai
    BGZF *fp = bgzf_open(args->ref_fai_fname,"r");
    if ( !fp ) error("Failed to read: %s\n", args->ref_fai_fname);
    kstring_t str = {0,0,0};
    while ( bgzf_getline(fp, '\n', &str) > 0 )
    {
        char keep, *tmp, *ptr = str.s;
        while ( *ptr && !isspace(*ptr) ) ptr++;
        if ( !*ptr ) error("Could not parse %s: %s\n", args->ref_fai_fname,str.s);
        keep = *ptr; *ptr = 0;
        if ( khash_str2int_has_key(regs->chr2i,str.s) ) continue;

        regs->nchr++;
        regs->chr = (chr_t*) realloc(regs->chr,sizeof(*regs->chr)*regs->nchr);
        chr_t *chr = &regs->chr[regs->nchr-1];
        memset(chr, 0, sizeof(chr_t));
        chr->name = strdup(str.s);
        *ptr = keep;
        chr->len  = strtol(ptr, &tmp, 10);
        if ( tmp==ptr ) error("Could not parse %s: %s\n", args->ref_fai_fname,str.s);
        khash_str2int_set(regs->chr2i, strdup(chr->name), regs->nchr-1);
        genome_len += chr->len;
    }
    if ( bgzf_close(fp)!=0 ) error("close failed: %s\n", args->ref_fai_fname);

    // read background+target regions, merge overlaps, remove duplicates, create a regidx
    if ( args->accessible_fname )
    {
        read_regions(args->accessible_fname,regs);
        read_regions(args->labeled_fname,regs);
        merge_regions(regs);
        args->bg_idx = regidx_init(NULL,NULL,NULL,0,NULL);
        for (i=0; i<regs->nreg; i++)
        {
            reg_t *reg = &regs->reg[i];
            char *chr_beg = regs->chr[reg->ichr].name;
            char *chr_end = chr_beg + strlen(chr_beg) - 1;
            regidx_push(args->bg_idx, chr_beg, chr_end, reg->beg, reg->end, NULL);
            if ( args->debug )
                fprintf(args->out_fh, "BG\t%s\t%d\t%d\n",chr_beg,reg->beg+1,reg->end+1);
            bg_len += reg->end - reg->beg + 1;
        }
        regs->nreg = 0;
        trim_chrs(regs,args->bg_idx);
    }

    // read the target regions, merge overlaps, remove duplicates, create a regidx
    regs->lbl2i = khash_str2int_init();
    read_regions(args->labeled_fname,regs);
    merge_regions(regs);
    args->tgt_idx = regidx_init(NULL,NULL,NULL,sizeof(int),NULL);
    for (i=0; i<regs->nreg; i++)
    {
        reg_t *reg = &regs->reg[i];
        char *chr_beg = regs->chr[reg->ichr].name;
        char *chr_end = chr_beg + strlen(chr_beg) - 1;
        regidx_push(args->tgt_idx, chr_beg, chr_end, reg->beg, reg->end, (void*)&reg->ilbl);
        if ( args->debug )
            fprintf(args->out_fh, "TGT\t%s\t%d\t%d\n",chr_beg,reg->beg+1,reg->end+1);
        if ( !args->bg_idx ) bg_len += reg->end - reg->beg + 1;
    }
    if ( !args->bg_idx ) trim_chrs(regs,args->tgt_idx);
    args->tgt_itr = regitr_init(args->tgt_idx);
    free(regs->reg);

    // initialize tracking arrays
    args->nobs = (uint32_t*) calloc(regs->nlbl,sizeof(*args->nobs));
    args->neq  = (uint32_t*) calloc(regs->nlbl,sizeof(*args->nexc));
    args->nexc = (uint32_t*) calloc(regs->nlbl,sizeof(*args->nexc));
    args->nfew = (uint32_t*) calloc(regs->nlbl,sizeof(*args->nfew));
    args->nhit = (uint32_t*) calloc(regs->nlbl,sizeof(*args->nhit));
    args->hit  = kbs_init(regs->nlbl);

    // read the calls
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
        int has_overlap = regidx_overlap(args->bg_idx ? args->bg_idx : args->tgt_idx, chr_beg,beg,end, NULL);
        if ( !has_overlap || end - beg + 1 > args->max_call_len )
        {
            if ( args->debug )
                fprintf(args->out_fh, "SKIP\t%s\t%d\t%d\n",chr_beg,beg+1,end+1);
            continue;
        }
        if ( args->debug )
            fprintf(args->out_fh, "CALL\t%s\t%d\t%d\n",chr_beg,beg+1,end+1);

        // which genes were hit? this is to initalize nobs
        kbs_clear(args->hit);
        if ( regidx_overlap(args->tgt_idx,chr_beg,beg,end,args->tgt_itr) )
        {
            // mark all genes hit by the call
            while ( regitr_overlap(args->tgt_itr) )
                kbs_insert(args->hit,regitr_payload(args->tgt_itr,int));

            kbitset_iter_t itr;
            kbs_start(&itr);
            while ((i = kbs_next(args->hit, &itr)) >= 0) args->nobs[i]++;
        }

        args->ncall++;
        args->call = (uint32_t*) realloc(args->call, sizeof(*args->call)*args->ncall);
        if ( !args->call ) error("Could not alloc %zu bytes\n", sizeof(*args->call)*args->ncall);
        args->call[args->ncall-1] = end - beg + 1;
    }
    if ( bgzf_close(fp)!=0 ) error("close failed: %s\n", args->calls_fname);
    free(str.s);

    if ( !args->ncall )
        error("Error: none of the calls intersects the tgt or bg regions (calls larger than `-m %d` were excluded)\n", args->max_call_len);

    args->dist = (uint32_t*) calloc(DIST_N(args),sizeof(*args->nhit));

    // sanity check for random placements: expected maximum number of attempts in inaccessible regions
    args->ntry = 10 * genome_len / bg_len;
}
static void destroy_data(args_t *args)
{
    regs_t *regs = &args->regs;
    int i;
    for (i=0; i<regs->nchr; i++) free(regs->chr[i].name);
    free(regs->chr);
    khash_str2int_destroy_free(regs->chr2i);
    for (i=0; i<regs->nlbl; i++) free(regs->lbl[i]);
    free(regs->lbl);
    khash_str2int_destroy(regs->lbl2i);
    if ( args->bg_idx ) regidx_destroy(args->bg_idx);
    regidx_destroy(args->tgt_idx);
    regitr_destroy(args->tgt_itr);
    free(args->nobs);
    free(args->neq);
    free(args->nexc);
    free(args->nfew);
    free(args->nhit);
    free(args->call);
    free(args->dist);
    kbs_destroy(args->hit);
}

static void run_test(args_t *args)
{
    regs_t *regs = &args->regs;
    memset(args->nhit,0,sizeof(*args->nhit)*regs->nlbl);

    int i, j;
    for (i=0; i<args->ncall; i++)
    {
        uint32_t len = args->call[i];
        kbs_clear(args->hit);
        int is_hit = 0, ntry;
        for (ntry=0; ntry<args->ntry; ntry++)
        {
            // randomly assign a chromosome
            int ichr = (double)random()/((double)RAND_MAX+1) * regs->nchr;
            chr_t *chr = &regs->chr[ichr];
            if ( chr->len < len ) continue;
            uint32_t pos = (double)random()/((double)RAND_MAX+1) * (chr->len - len + 1);
            is_hit = regidx_overlap(args->tgt_idx,chr->name,pos,pos+len-1,args->tgt_itr);
            if ( is_hit ) break;
            if ( args->bg_idx && regidx_overlap(args->bg_idx,chr->name,pos,pos+len-1,NULL) ) break;
        }
        if ( ntry==args->ntry ) error("Check me: this does not look right, too many failures\n");
        if ( !is_hit ) continue;

        // mark all genes hit by the call
        while ( regitr_overlap(args->tgt_itr) )
            kbs_insert(args->hit,regitr_payload(args->tgt_itr,int));
        
        kbitset_iter_t itr;
        kbs_start(&itr);
        while ((j = kbs_next(args->hit, &itr)) >= 0) args->nhit[j]++;
    }

    // evaluate if a gene was hit more/fewer times than observed
    for (i=0; i<regs->nlbl; i++)
    {
        if ( args->nhit[i] == args->nobs[i] ) args->neq[i]++;
        else if ( args->nhit[i] > args->nobs[i] ) args->nexc[i]++;
        else if ( args->nhit[i] < args->nobs[i] ) args->nfew[i]++;

        int idx = DIST_IDX(args,i,args->nhit[i]);
        args->dist[idx]++;
    }
}

static void usage(void)
{
    error(
        "\n"
        "Program: recurrence-test %s\n"
        "License: The MIT/Expat license\n"
        "This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        "About: Run recurrence enrichment (permutation) test.\n"
        "Usage: recurrence-test [OPTIONS]\n"
        "Options:\n"
        "   -a, --accessible-regs FILE      Optional list of accessible regions (otherwise equals to -l): chr,beg,end\n"
        "   -c, --calls FILE                Calls: chr,beg,end\n"
        "   -d, --debug-regions             Print the spliced regions (and stop)\n"
        "   -f, --ref-fai FILE              Chromosome lengths, given for example as fai index: chr,length\n"
        "   -l, --labeled-regs FILE         Labeled regions to test the recurrence: chr,beg,end,label\n"
        "   -m, --max-call-length INT       Skip big calls [10e6]\n"
        "   -n, --niter INT                 Number of iterations to run\n"
        "   -o, --output FILE               Place output in FILE\n"
        "   -s, --random-seed INT           Random seed\n"
        "\n"
        "Examples:\n"
        "   [todo]"
        "\n",
        UTILS_VERSION
        );
}

int main(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->niter = 0;
    args->max_call_len = 10e6;
    static struct option loptions[] =
    {
        {"max-call-len",required_argument,NULL,'m'},
        {"ref-fai",required_argument,NULL,'f'},
        {"output",required_argument,NULL,'o'},
        {"ninter",required_argument,NULL,'n'},
        {"debug-regions",required_argument,NULL,'d'},
        {"calls",required_argument,NULL,'c'},
        {"accessible-regs",required_argument,NULL,'b'},
        {"labeled-regs",required_argument,NULL,'l'},
        {"random-seed",required_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int i,j, c, seed = -1;
    while ((c = getopt_long(argc, argv, "?hc:a:l:n:ds:f:o:m:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'm': 
                args->max_call_len = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --max-call-length %s\n", optarg);
                break;
            case 's': 
                seed = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --random-seed %s\n", optarg);
                break;
            case 'n': 
                args->niter = strtod(optarg, &tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse --niter%s\n", optarg);
                break;
            case 'f': args->ref_fai_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'c': args->calls_fname = optarg; break;
            case 'l': args->labeled_fname = optarg; break;
            case 'a': args->accessible_fname = optarg; break;
            case 'd': args->debug++; break;
            case 'h':
            case '?':
            default: usage(); break;
        }
    }
    if ( argc < 2 ) usage();
    if ( !args->ref_fai_fname ) error("Missing the -f option\n");
    if ( !args->calls_fname ) error("Missing the -c option\n");
    if ( !args->labeled_fname ) error("Missing the -l option\n");

    args->out_fh = args->output_fname ? fopen(args->output_fname,"w") : stdout;
    if ( !args->out_fh ) error("Error: Could not open %s for writing\n", args->output_fname);
    if ( seed==-1 )
    {
        struct timeval time; 
        gettimeofday(&time,NULL);
        seed = time.tv_sec + time.tv_usec;
    }
    srandom(seed);

    init_data(args);

    fprintf(args->out_fh, "VERSION\t%s\n",UTILS_VERSION);
    fprintf(args->out_fh, "CMD\t%s",argv[0]); for (c=1; c<argc; c++) fprintf(args->out_fh, " %s", argv[c]); fprintf(args->out_fh, "\n");
    fprintf(args->out_fh, "SEED\t%d\n", seed);
    fprintf(args->out_fh, "NITER\t%e\n", (double)args->niter);
    fprintf(args->out_fh,"# TEST:\n");
    fprintf(args->out_fh,"#     - gene .. the region label\n");
    fprintf(args->out_fh,"#     - nObserved .. number of hits observed in input data in this gene\n");
    fprintf(args->out_fh,"#     - nSimFewer .. number of simulations with fewer hits in the gene than observed in the data\n");
    fprintf(args->out_fh,"#     - nSimEqual .. as above, but with the same number of hits than observed\n");
    fprintf(args->out_fh,"#     - nSimExceeded .. as above, but with the more hits than observed\n");
    fprintf(args->out_fh,"# DIST:\n");
    fprintf(args->out_fh,"#     - gene\n");
    fprintf(args->out_fh,"#     - nSimHits .. the number of simulations with 0,1,2,etc. hits in the gene\n");

    for (i=0; i<args->niter; i++) run_test(args);

    regs_t *regs = &args->regs;
    fprintf(args->out_fh,"# [1]TEST\t[2]gene\t[3]nObserved\t[4]nSimFewer\t[5]nSimEqual\t[6]nSimExceeded\n");
    for (i=0; i<regs->nlbl; i++)
        fprintf(args->out_fh,"TEST\t%s\t%u\t%u\t%u\t%u\n",regs->lbl[i],args->nobs[i],args->nfew[i],args->neq[i],args->nexc[i]);

    fprintf(args->out_fh,"# [1]DIST\t[2]gene\t[3-]nSimHits (0,1,etc. hits)\n");
    for (i=0; i<regs->nlbl; i++)
    {
        fprintf(args->out_fh,"DIST\t%s",regs->lbl[i]);
        for (j=0; j<=args->ncall; j++)
        {
            int idx = DIST_IDX(args,i,j);
            fprintf(args->out_fh,"\t%u",args->dist[idx]);
        }
        fprintf(args->out_fh,"\n");
    }

    destroy_data(args);
    free(args);

    return 0;
}

