/* 
    Copyright (C) 2018-2019 Genome Research Ltd.
    
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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/khash_str2int.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <getopt.h>
#include <assert.h>
#include "regidx.h"
#include "version.h"

#define ANN_NBP     1
#define ANN_FRAC    2

typedef struct
{
    int n,m;
    char **off, *rmme;
}
cols_t;

typedef struct
{
    void *name2idx;
    cols_t *cols, *annots;
    int dummy;
}
hdr_t;

typedef struct
{
    char *fname;
    hdr_t *hdr;
    cols_t *core, *match, *transfer, *annots;
    int *core_idx, *match_idx, *transfer_idx, *annots_idx;
    int grow_n;
}
dat_t;

// This is for the special -a annotations, keeps a list of 
// soure regions that hit the destination region. The start
// coordinates are converted to beg<<1 and end coordinates
// to (end<<1)+1.
#define NBP_SET_BEG(x) ((x)<<1)
#define NBP_SET_END(x) (((x)<<1)+1)
#define NBP_GET(x)     ((x)>>1)
#define NBP_IS_BEG(x)  (((x)&1)==0)
#define NBP_IS_END(x)  (((x)&1)==1)
typedef struct
{
    int n,m;
    uint32_t *regs; // change to uint64_t for very large genomes
}
nbp_t;

typedef struct
{
    nbp_t *nbp;
    dat_t dst, src;
    char *core_str, *match_str, *transfer_str, *annots_str;
    char *temp_dir;
    int allow_dups, reciprocal;
    float overlap;
    regidx_t *idx;
    regitr_t *itr;
    htsFile *fh;
    kstring_t tmp_kstr;
    cols_t *tmp_cols;               // the -t transfer fields to write for each line
    khash_t(str2int) **tmp_hash;    // lookup tables for tmp_cols
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

static nbp_t *nbp_init(void)
{
    return calloc(1,sizeof(nbp_t));
}
static void nbp_destroy(nbp_t *nbp)
{
    free(nbp->regs);
    free(nbp);
}
static inline void nbp_reset(nbp_t *nbp)
{
    nbp->n = 0;
}
static inline void nbp_add(nbp_t *nbp, uint32_t beg, uint32_t end)
{
    if ( end >= REGIDX_MAX>>1 ) error("Error: the coordinate is too big (%u). Possible todo: switch to uint64_t\n",end);
    nbp->n += 2;
    if ( nbp->n >= nbp->m )
    {
        nbp->m += 2;
        nbp->regs = realloc(nbp->regs, nbp->m*sizeof(*nbp->regs));
    }
    nbp->regs[nbp->n - 2] = NBP_SET_BEG(beg);
    nbp->regs[nbp->n - 1] = NBP_SET_END(end);
}
static int compare_uint32(const void *aptr, const void *bptr)
{
    uint32_t a = *(const uint32_t*) aptr;
    uint32_t b = *(const uint32_t*) bptr;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}
static int nbp_length(nbp_t *nbp)
{
    qsort(nbp->regs, nbp->n, sizeof(*nbp->regs), compare_uint32);
    int i, nopen = 0, length = 0;
    uint32_t beg = 0;
    for (i=0; i<nbp->n; i++)
    {
        if ( NBP_IS_BEG(nbp->regs[i]) )
        {
            if ( !nopen ) beg = NBP_GET(nbp->regs[i]);
            nopen++;
        }
        else nopen--;
        assert( nopen>=0 );
        if ( nopen==0 && beg>0 ) length += NBP_GET(nbp->regs[i]) - beg + 1;
    }
    return length;
}

cols_t *cols_split(const char *line, cols_t *cols, char delim)
{
    if ( !cols ) cols = (cols_t*) calloc(1,sizeof(cols_t));
    if ( cols->rmme ) free(cols->rmme);
    cols->n = 0;
    cols->rmme = strdup(line);
    char *ss = cols->rmme;
    while (1)
    {
        char *se = ss;
        while ( *se && *se!=delim ) se++;
        char tmp = *se;
        *se = 0;
        cols->n++;
        if ( cols->n > cols->m )
        {
            cols->m += 10;
            cols->off = realloc(cols->off, sizeof(*cols->off)*cols->m);
        }
        cols->off[ cols->n - 1 ] = ss;
        if ( !tmp ) break;
        ss = se + 1;
    }
    return cols;
}
// Can be combined with cols_split() but is much slower.
// The string must exist throughout the life of cols unless initialized with cols_split().
void cols_append(cols_t *cols, char *str)
{
    if ( cols->rmme )
    {
        size_t str_len = strlen(str);
        size_t lst_len = strlen(cols->off[ cols->n - 1 ]);
        size_t tot_len = 2 + str_len + lst_len + (cols->off[ cols->n - 1 ] - cols->rmme);

        cols_t *tmp_cols = (cols_t*)calloc(1,sizeof(cols_t));
        tmp_cols->rmme = calloc(tot_len,1);
        tmp_cols->off  = calloc(cols->n+1,sizeof(*tmp_cols->off));

        char *ptr = tmp_cols->rmme;
        int i;
        for (i=0; i<cols->n; i++)
        {
            size_t len = strlen(cols->off[i]);
            memcpy(ptr, cols->off[i], len);
            tmp_cols->off[i] = ptr;
            ptr += len + 1;
        }
        memcpy(ptr, str, str_len);
        tmp_cols->off[i] = ptr;

        free(cols->off);
        free(cols->rmme);
        cols->rmme = tmp_cols->rmme;
        cols->off  = tmp_cols->off;
        cols->n    = cols->n+1;
        cols->m    = cols->n;
        free(tmp_cols);
        return;
    }
    cols->n++;
    if ( cols->n > cols->m )
    {
        cols->m++;
        cols->off = realloc(cols->off,sizeof(*cols->off)*cols->m);
    }
    cols->off[cols->n-1] = str;
}
void cols_clear(cols_t *cols)
{
    if ( !cols ) return;
    free(cols->rmme);
    free(cols->off);
    cols->rmme = NULL;
    cols->off  = NULL;
}
void cols_destroy(cols_t *cols)
{
    if ( !cols ) return;
    cols_clear(cols);
    free(cols);
}

int parse_tab_with_payload(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    if ( line[0]=='#' )
    {
        *((cols_t**)payload) = NULL;
        return -1;
    }

    dat_t *dat = (dat_t*) usr;

    cols_t *cols = cols_split(line, NULL, '\t');
    *((cols_t**)payload) = cols;

    if ( cols->n <= dat->core_idx[0] ) error("Could not parse %s\n", line);
    *chr_beg = cols->off[ dat->core_idx[0] ];
    *chr_end = *chr_beg + strlen(*chr_beg) - 1;

    if ( cols->n <= dat->core_idx[1] ) error("Could not parse %s\n", line);
    char *tmp, *ptr = cols->off[ dat->core_idx[1] ];
    *beg = strtod(ptr, &tmp);
    if ( tmp==ptr ) error("Error: Could not parse: %s\n", line);

    if ( cols->n <= dat->core_idx[2] ) error("Could not parse %s\n", line);
    ptr = cols->off[ dat->core_idx[2] ];
    *end = strtod(ptr, &tmp);
    if ( tmp==ptr ) error("Error: Could not parse: %s\n", line);

    return 0;
}
void free_payload(void *payload)
{
    cols_t *cols = *((cols_t**)payload);
    cols_destroy(cols);
}

hdr_t *parse_header(char *fname)
{
    hdr_t *hdr = (hdr_t*) calloc(1, sizeof(hdr_t));

    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) error("Failed to open: %s\n", fname);
    kstring_t line = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &line) <= 0 ) error("Failed to read: %s\n", fname);

    cols_t *cols = cols_split(line.s, NULL, '\t');
    assert(cols && cols->n);
    if ( cols->off[0][0] != '#' )
    {
        kstring_t str = {0,0,0};
        int i, n = cols->n;
        for (i=0; i<n; i++)
        {
            if ( i>0 ) kputc('\t', &str);
            kputw(i+1, &str);
        }
        cols_destroy(cols);
        cols = cols_split(str.s, NULL, '\t');
        free(str.s);
        hdr->dummy = 1;
    }

    hdr->name2idx = khash_str2int_init();
    int i;
    for (i=0; i<cols->n; i++)
    {
        char *ss = cols->off[i];
        while ( *ss && (*ss=='#' || isspace(*ss)) ) ss++;
        if ( !*ss ) error("Could not parse the header field \"%s\": %s\n", cols->off[i],line.s);
        if ( *ss=='[' )
        {
            char *se = ss+1;
            while ( *se && isdigit(*se) ) se++;
            if ( *se==']' ) ss = se + 1;
        }
        while ( *ss && (*ss=='#' || isspace(*ss)) ) ss++;
        if ( !*ss ) error("Could not parse the header field \"%s\": %s\n", cols->off[i],line.s);
        cols->off[i] = ss;
        khash_str2int_set(hdr->name2idx, cols->off[i], i);
    }
    hdr->cols = cols;

    free(line.s);
    hts_close(fp);

    return hdr;
}
void write_header(hdr_t *hdr)
{
    if ( hdr->dummy ) return;
    int i;
    kstring_t str = {0,0,0};
    kputc('#', &str);
    for (i=0; i<hdr->cols->n; i++)
    {
        if ( i>0 ) kputc('\t', &str);
        ksprintf(&str,"[%d]", i+1);
        kputs(hdr->cols->off[i], &str);
    }
    if ( hdr->annots )
    {
        for (i=0; i<hdr->annots->n; i++)
        {
            if ( str.l > 1 ) kputc('\t', &str);
            kputs(hdr->annots->off[i], &str);
        }
    }
    kputc('\n',&str);
    if ( fwrite(str.s, str.l, 1, stdout) != 1 ) error("Failed to write %d bytes\n", str.l);
    free(str.s);
}
void destroy_header(hdr_t *hdr)
{
    if ( hdr->cols ) cols_destroy(hdr->cols);
    khash_str2int_destroy(hdr->name2idx);
    free(hdr);
}

void sanity_check_columns(char *fname, hdr_t *hdr, cols_t *cols, int **col2idx, int force)
{
    *col2idx = (int*)malloc(sizeof(int)*cols->n);
    int i, idx;
    for (i=0; i<cols->n; i++)
    {
        if ( khash_str2int_get(hdr->name2idx, cols->off[i], &idx) < 0 )
        {
            if ( !force ) error("The key \"%s\" not found in %s\n", cols->off[i],fname);
            idx = -1;
        }
        (*col2idx)[i] = idx;
    }
}
void init_data(args_t *args)
{
    args->dst.hdr = parse_header(args->dst.fname);
    args->src.hdr = parse_header(args->src.fname);

    // -c, core columns
    if ( !args->core_str ) args->core_str = "chr,beg,end:chr,beg,end"; 
    cols_t *tmp = cols_split(args->core_str, NULL, ':');
    args->src.core = cols_split(tmp->off[0],NULL,',');
    args->dst.core = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
    sanity_check_columns(args->src.fname, args->src.hdr, args->src.core, &args->src.core_idx, 0);
    sanity_check_columns(args->dst.fname, args->dst.hdr, args->dst.core, &args->dst.core_idx, 0);
    if ( args->src.core->n!=3 || args->dst.core->n!=3 ) error("Expected three columns: %s\n", args->core_str);

    // -m, match columns
    if ( args->match_str )
    {
        tmp = cols_split(args->match_str, tmp, ':');
        args->src.match = cols_split(tmp->off[0],NULL,',');
        args->dst.match = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
        sanity_check_columns(args->src.fname, args->src.hdr, args->src.match, &args->src.match_idx, 0);
        sanity_check_columns(args->dst.fname, args->dst.hdr, args->dst.match, &args->dst.match_idx, 0);
        if ( args->src.match->n != args->dst.match->n ) error("Expected equal number of columns: %s\n", args->match_str);
    }

    // -t, transfer columns
    if ( !args->transfer_str ) error("Missing the -t, --transfer option\n");
    tmp = cols_split(args->transfer_str, tmp, ':');
    args->src.transfer = cols_split(tmp->off[0],NULL,',');
    args->dst.transfer = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
    sanity_check_columns(args->src.fname, args->src.hdr, args->src.transfer, &args->src.transfer_idx, 1);
    sanity_check_columns(args->dst.fname, args->dst.hdr, args->dst.transfer, &args->dst.transfer_idx, 1);
    if ( args->src.transfer->n != args->dst.transfer->n ) error("Expected equal number of columns: %s\n", args->transfer_str);
    int i;
    for (i=0; i<args->src.transfer->n; i++)
    {
        if ( args->src.transfer_idx[i]==-1 )
        {
            cols_append(args->src.hdr->cols,args->src.transfer->off[i]);
            args->src.transfer_idx[i] = -args->src.hdr->cols->n;    // negative index indicates different ptr location
            args->src.grow_n++;
        }
    }
    for (i=0; i<args->dst.transfer->n; i++)
    {
        if ( args->dst.transfer_idx[i]==-1 )
        {
            cols_append(args->dst.hdr->cols,args->dst.transfer->off[i]);
            args->dst.transfer_idx[i] = args->dst.hdr->cols->n - 1;
            args->dst.grow_n++;
        }
    }
    args->tmp_cols = (cols_t*)calloc(args->src.transfer->n,sizeof(cols_t));
    args->tmp_hash = (khash_t(str2int)**)calloc(args->src.transfer->n,sizeof(khash_t(str2int)*));
    for (i=0; i<args->src.transfer->n; i++)
        args->tmp_hash[i] = khash_str2int_init();
    cols_destroy(tmp);

    // -a, annotation columns
    if ( args->annots_str )
    {
        tmp = cols_split(args->annots_str, NULL, ':');
        args->src.annots = cols_split(tmp->off[0],NULL,',');
        args->dst.annots = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
        if ( args->src.annots->n!=args->dst.annots->n ) error("Different number of src and dst columns in %s\n",args->annots_str);
        args->dst.annots_idx = (int*) malloc(sizeof(int)*args->dst.annots->n);
        for (i=0; i<args->src.annots->n; i++)
        {
            if ( !strcasecmp(args->src.annots->off[i],"nbp") ) args->dst.annots_idx[i] = ANN_NBP;
            else if ( !strcasecmp(args->src.annots->off[i],"frac") ) args->dst.annots_idx[i] = ANN_FRAC;
            else error("The annotation \"%s\" is not recognised\n", args->src.annots->off[i]);
        }
        args->nbp = nbp_init();
    }

    args->idx = regidx_init(args->src.fname, parse_tab_with_payload,free_payload,sizeof(cols_t),&args->src);
    args->itr = regitr_init(args->idx);

    args->fh = hts_open(args->dst.fname,"r");
    if ( !args->fh ) error("Failed to read: %s\n", args->dst.fname);
}
void destroy_data(args_t *args)
{
    if ( hts_close(args->fh)!=0 ) error("Failed to close: %s\n", args->dst.fname);
    int i;
    for (i=0; i<args->src.transfer->n; i++)
        khash_str2int_destroy(args->tmp_hash[i]);
    free(args->tmp_hash);
    for (i=0; i<args->src.transfer->n; i++) cols_clear(&args->tmp_cols[i]);
    free(args->tmp_cols);
    cols_destroy(args->src.core);
    cols_destroy(args->dst.core);
    cols_destroy(args->src.match);
    cols_destroy(args->dst.match);
    cols_destroy(args->src.transfer);
    cols_destroy(args->dst.transfer);
    if ( args->src.annots ) cols_destroy(args->src.annots);
    if ( args->dst.annots ) cols_destroy(args->dst.annots);
    if ( args->nbp ) nbp_destroy(args->nbp);
    destroy_header(args->src.hdr);
    destroy_header(args->dst.hdr);
    free(args->src.core_idx);
    free(args->dst.core_idx);
    free(args->src.match_idx);
    free(args->dst.match_idx);
    free(args->src.transfer_idx);
    free(args->dst.transfer_idx);
    free(args->src.annots_idx);
    free(args->dst.annots_idx);
    if (args->itr) regitr_destroy(args->itr);
    if (args->idx) regidx_destroy(args->idx);
    free(args->tmp_kstr.s);
}

static inline void write_string(args_t *args, char *str, size_t len)
{
    if ( len==0 ) len = strlen(str);
    if ( fwrite(str, len, 1, stdout) != 1 ) error("Failed to write %d bytes\n", len);
}

void process_line(args_t *args, char *line, size_t size)
{
    char *chr_beg, *chr_end;
    uint32_t beg, end;
    cols_t *dst_cols = NULL;
    int i,j;
    int ret = parse_tab_with_payload(line, &chr_beg, &chr_end, &beg, &end, &dst_cols, &args->dst);
    if ( ret==-1 )
    {
        cols_destroy(dst_cols);
        return;
    }

    if ( !regidx_overlap(args->idx, chr_beg,beg,end, args->itr) )
    {
        write_string(args, line, size);
        write_string(args, "\n", 1);
        cols_destroy(dst_cols);
        return;
    }

    for (i=0; i<args->src.transfer->n; i++)
    {
        args->tmp_cols[i].n = 0;
        kh_clear(str2int, args->tmp_hash[i]);
    }

    if ( args->nbp ) nbp_reset(args->nbp);

    int has_overlap = 0;
    while ( regitr_overlap(args->itr) )
    {
        if ( args->overlap )
        {
            float len1 = end - beg + 1;
            float len2 = args->itr->end - args->itr->beg + 1;
            float isec = (args->itr->end < end ? args->itr->end : end) - (args->itr->beg > beg ? args->itr->beg : beg) + 1;
            if ( args->reciprocal )
            {
                if ( isec/len1 < args->overlap || isec/len2 < args->overlap ) continue;
            }
            else
            {
                if ( isec/len1 < args->overlap && isec/len2 < args->overlap ) continue;
            }
        }
        cols_t *src_cols = regitr_payload(args->itr,cols_t*);
        if ( args->dst.match && args->dst.match->n )
        {
            for (i=0; i<args->dst.match->n; i++)
            {
                if ( args->dst.match_idx[i] > dst_cols->n ) error("Could not parse: %s\n", line);
                char *dst = dst_cols->off[ args->dst.match_idx[i] ];
                char *src = src_cols->off[ args->src.match_idx[i] ];
                if ( strcmp(dst,src) ) break;
            }
            if ( i != args->dst.match->n ) continue;
        }

        if ( args->nbp )
            nbp_add(args->nbp, args->itr->beg >= beg ? args->itr->beg : beg, args->itr->end <= end ? args->itr->end : end);

        for (i=0; i<args->src.transfer->n; i++)
        {
            char *str;
            if ( args->src.transfer_idx[i] >= 0 )
                str = src_cols->off[ args->src.transfer_idx[i] ];                   // transfer a value from the src file
            else
                str = args->src.hdr->cols->off[ -args->src.transfer_idx[i] - 1 ];    // non-existent field in src, use a default value

            if ( !args->allow_dups )
            {
                if ( khash_str2int_has_key(args->tmp_hash[i],str) ) continue;
                khash_str2int_set(args->tmp_hash[i],str,1);
            }
            cols_append(&args->tmp_cols[i], str);
            has_overlap += strlen(str);
        }
    }

    if ( !has_overlap )
    {
        write_string(args, line, size);
        write_string(args, "\n", 1);
        cols_destroy(dst_cols);
        return;
    }

    size_t len;
    args->tmp_kstr.l = 0;
    ks_resize(&args->tmp_kstr, has_overlap*3 + args->src.transfer->n*2);
    for (i=0; i<args->src.transfer->n; i++)
    {
        char *off = dst_cols->off[ args->dst.transfer_idx[i] ] = args->tmp_kstr.s + args->tmp_kstr.l;
        cols_t *ann = &args->tmp_cols[i];
        if ( !ann->n ) { off[0] = '.'; off[1] = 0; args->tmp_kstr.l += 2; continue; }
        for (j=0; j<ann->n; j++)
        {
            if ( j>0 ) { off[0] = ','; off++; args->tmp_kstr.l++; }
            len = strlen(ann->off[j]);
            memcpy(off, ann->off[j], len);
            off += len;
            args->tmp_kstr.l += len;
        }
        off[0] = 0;
        args->tmp_kstr.l++;
    }
    write_string(args, dst_cols->off[0], 0);
    for (i=1; i<dst_cols->n; i++)
    {
        write_string(args, "\t", 1);
        write_string(args, dst_cols->off[i], 0);
    }

    if ( args->src.annots )
    {
        args->tmp_kstr.l = 0;
        int len = nbp_length(args->nbp);
        for (i=0; i<args->dst.annots->n; i++)
        {
            if ( args->dst.annots_idx[i]==ANN_NBP ) 
            {
                kputc('\t',&args->tmp_kstr);
                kputw(len,&args->tmp_kstr);
            }
            else if ( args->dst.annots_idx[i]==ANN_FRAC ) 
            {
                kputc('\t',&args->tmp_kstr);
                kputd((double)len/(end-beg+1),&args->tmp_kstr);
            }
        }
        write_string(args, args->tmp_kstr.s, args->tmp_kstr.l);
    }

    write_string(args, "\n", 1);
    cols_destroy(dst_cols);
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Annotate regions in DST file with texts from overlapping regions in SRC file.\n"
        "       The transfer of annotations can be conditioned on matching values in one or more\n"
        "       columns (-m), multiple columns can be transferred (-t).\n"
        "       All indexes and coordinates are 1-based and inclusive.\n"
        "Usage: annot-regs [OPTIONS] DST\n"
        "Options:\n"
        "       --allow-dups                Add annotations multiple times\n"
        "   -a, --annotate list             Add special annotations:\n"
        "                                       frac .. fraction of the destination region with an overlap\n"
        "                                       nbp  .. number of source base pairs in the overlap\n"
        "   -c, --core src:dst              Core columns [chr,beg,end:chr,beg,end]\n"
        "   -d, --dst-file file             Destination file\n"
        "   -m, --match src:dst             Require match in these columns\n"
        "   -o, --overlap float             Minimum required overlap (non-reciprocal, unless -r is given)\n"
        "   -r, --reciprocal                Require reciprocal overlap\n"
        "   -s, --src-file file             Source file\n"
        "   -t, --transfer src:dst          Columns to transfer. If src column does not exist, interpret\n"
        "                                   as the default value to use. If the dst column does not exist,\n"
        "                                   a new column is created.\n"
        "       --version                   Print version string and exit\n"
        "Examples:\n"
        "   # Header is present, match and transfer by column name\n"
        "   annot-regs -s src.txt.gz -d dst.txt.gz -c chr,beg,end:chr,beg,end -m type,sample:type,smpl -t tp/fp:tp/fp\n"
        "\n"
        "   # Header is not present, match and transfer by column index (1-based)\n"
        "   annot-regs -s src.txt.gz -d dst.txt.gz -c 1,2,3:1,2,3 -m 4,5:4,5 -t 6:6\n"
        "\n"
        "   # If the dst part is not given, the program assumes that the src:dst columns are identical\n"
        "   annot-regs -s src.txt.gz -d dst.txt.gz -c chr,beg,end -m type,sample -t tp/fp\n"
        "\n";
}

int main(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    static struct option loptions[] =
    {
        {"allow-dups",no_argument,NULL,0},
        {"version",no_argument,NULL,1},
        {"annotate",required_argument,NULL,'a'},
        {"core",required_argument,NULL,'c'},
        {"dst-file",required_argument,NULL,'d'},
        {"match",required_argument,NULL,'m'},
        {"overlap",required_argument,NULL,'o'},
        {"reciprocal",no_argument,NULL,'r'},
        {"src-file",required_argument,NULL,'s'},
        {"transfer",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int c;
    while ((c = getopt_long(argc, argv, "hc:d:m:o:s:t:T:ra:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  0 : args->allow_dups = 1; break;
            case  1 : printf("%s\n",ANNOT_REGS_VERSION); return 0; break;
            case 'r': args->reciprocal = 1; break;
            case 'c': args->core_str  = optarg; break;
            case 'd': args->dst.fname = optarg; break;
            case 'm': args->match_str = optarg; break;
            case 'a': args->annots_str = optarg; break;
            case 'o': 
                args->overlap = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --overlap %s\n", optarg);
                break;
            case 's': args->src.fname = optarg; break;
            case 't': args->transfer_str = optarg; break;
            case 'h':
            case '?':
            default: error("%sVersion: %s\n\n", usage_text(), ANNOT_REGS_VERSION); break;
        }
    }
    if ( argc==1 ) error("%sVersion: %s\n\n",usage_text(), ANNOT_REGS_VERSION);
    if ( !args->dst.fname ) error("Missing the -d option\n");
    if ( !args->src.fname ) error("Missing the -s option\n");

    init_data(args);
    kstring_t line = {0,0,0};

    write_header(args->dst.hdr);
    while ( hts_getline(args->fh, KS_SEP_LINE, &line) > 0 )
    {
        int i;
        for (i=0; i<args->dst.grow_n; i++) kputs("\t.", &line);
        process_line(args, line.s, line.l);
    }

    free(line.s);
    destroy_data(args);
    free(args);

    return 0;
}

