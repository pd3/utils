/* 
    Copyright (C) 2018 Genome Research Ltd.
    
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
#include "regidx.h"

typedef struct
{
    float min_overlap;
    int only_unique, new_column, require_match;
    char *fname, *src_annots, *default_miss, *default_hit;
    regidx_t *idx;
    regitr_t *itr;
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


// parse the DST file
int parse_tab(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, char **payload)
{
    *payload = NULL;

    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se )
    {
        // just the chromosome name
        *beg = 0;
        *end = REGIDX_MAX;
        return 0;
    }

    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) error("Error: Could not parse: %s\n", line);

    ss = se+1;
    *end = strtod(ss, &se);
    if ( ss==se ) error("Error: Could not parse: %s\n", line);

    if ( *beg > *end ) error("Error: beg>end: %s\n", line);

    if ( !*se ) return 0;

    // skip spaces (there should be no more than one)
    ss = se;
    while ( *ss && isspace(*ss) ) ss++;

    // read the payload
    se = *ss ? se + 1 : ss;
    while ( *se && !isspace(*se) ) se++;

    char tmp = *se; *se = 0;
    *payload = strdup(ss);
    *se = tmp;

    return 0;
}

// parse the SRC file with annotations, expect three or four columns
int regidx_parse_tab_with_payload(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    args_t *args = (args_t*) usr;

    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se )
    {
        // just the chromosome name
        *beg = 0;
        *end = REGIDX_MAX;
        return -1;
    }

    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) error("Error: Could not parse: %s\n", line);

    ss = se+1;
    *end = strtod(ss, &se);
    if ( ss==se ) error("Error: Could not parse: %s\n", line);

    if ( *beg > *end ) error("Error: beg>end: %s\n", line);

    // skip spaces (there should be no more than one)
    ss = se;
    while ( *ss && isspace(*ss) ) ss++;

    // read the payload
    se = *ss ? se + 1 : ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*ss )
    {
        if ( !args->default_hit ) error("Error: expected 4th colum or the -a option: %s\n", line);
        *((char**)payload) = strdup(args->default_hit);
        return 0;
    }
    *se = 0;
    *((char**)payload) = strdup(ss);

    return 0;
}
void free_payload(void *payload)
{
    free(*((char**)payload));
}


static const char *usage_text(void)
{
    return 
        "\n"
        "About: Annotate DST regions with texts from overlapping regions of SRC.\n"
        "       Three columns are expected in DST and SRC (chr,beg,end) and an optional\n"
        "       fourth column composed of comma-separated strings where matching annotations\n"
        "       are searched for. The DST file can have more columns, but only the fourth\n"
        "       is used for annotation matching.\n"
        "       Coordinates are 1-based and inclusive.\n"
        "Usage: annot-regs [OPTIONS] DST\n"
        "Options:\n"
        "   -a, --dflt-annots STR[,STR]     Default miss and hit annotation to add\n"
        "   -m, --min-overlap FLT           Minimum non-reciprocal overlap [any overlap]\n"
        "   -n, --new-column                Put annotations in a new column\n"
        "   -r, --require                   Require a SRC annotation to be already present in DST\n"
        "   -s, --source-annots SRC         Regions with source annotations\n"
        "   -u, --only-unique               Do not add the same annotation multiple times\n"
        "Example:\n"
        "   # Find overlapping regions appending the SRC annotation as a new column.\n"
        "   # If no overlap found, add a single dot.\n"
        "   annot-regs -a . -nu -s source.bed destination.bed\n"
        "\n"
        "   # Find overlapping regions with matching annotations, adding a new column with PASS\n"
        "   # if found and FAIL otherwise\n"
        "   annot-regs -a FAIL,PASS -nur -s source.bed destination.bed\n"
        "\n";
}

int main(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    static struct option loptions[] =
    {
        {"new-column",no_argument,NULL,'n'},
        {"dflt-annots",required_argument,NULL,'a'},
        {"only-unique",no_argument,NULL,'u'},
        {"require",no_argument,NULL,'r'},
        {"min-overlap",required_argument,NULL,'m'},
        {"source-annots",required_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int c, i;
    while ((c = getopt_long(argc, argv, "hs:ua:m:nr",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'r': args->require_match = 1; break;
            case 'm': 
                args->min_overlap = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --overlap %s\n", optarg);
                break;
            case 'u': args->only_unique = 1; break;
            case 'n': args->new_column = 1; break;
            case 's': args->src_annots = optarg; break;
            case 'a': 
                args->default_miss = optarg; 
                tmp = optarg;
                while ( *tmp && *tmp!=',' ) tmp++;
                if ( *tmp ) { *tmp = 0; args->default_hit = tmp+1; }
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( !args->src_annots ) error("%s",usage_text());
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    args->idx = regidx_init(args->src_annots,regidx_parse_tab_with_payload,free_payload,sizeof(char*),args);
    args->itr = regitr_init(args->idx);
    kstring_t str = {0,0,0}, line = {0,0,0};
    htsFile *fp = hts_open(args->fname,"r");
    if ( !fp ) error("Failed to read: %s\n", args->fname);

    while ( hts_getline(fp, KS_SEP_LINE, &line) > 0 )
    {
        char *chr_beg, *chr_end, *dst_payload;
        uint32_t beg, end;
        int ret = parse_tab(line.s, &chr_beg, &chr_end, &beg, &end, &dst_payload);
        if ( ret ) error("Failed to parse: %s\n", line.s);

        void *dst_annots_hash = NULL;

        // Parse existing annotations in DST, this is to prevent duplicate entries with multiple
        // overlapping records
        if ( args->only_unique )
        {
            // a short lived hash - in case of multiple overlapping regions, add each annotation only once
            dst_annots_hash = khash_str2int_init();

            if ( !args->new_column )
            {
                int ndst_annots = 0;
                char **dst_annots;
                dst_annots = hts_readlist(dst_payload, 0, &ndst_annots);
                for (i=0; i<ndst_annots; i++) khash_str2int_set(dst_annots_hash, strdup(dst_annots[i]), 1);
                for (i=0; i<ndst_annots; i++) free(dst_annots[i]);
                free(dst_annots);
            }
        }

        // Find the overlap
        str.l = 0;
        char tmp = chr_end[1]; chr_end[1] = 0;
        if ( regidx_overlap(args->idx, chr_beg,beg,end, args->itr) )
        {
            int first = 1;
            while ( regitr_overlap(args->itr) )
            {
                if ( args->min_overlap )
                {
                    float len1 = end - beg + 1;
                    float len2 = args->itr->end - args->itr->beg + 1;
                    float isec = (args->itr->end < end ? args->itr->end : end) - (args->itr->beg > beg ? args->itr->beg : beg) + 1;
                    if ( isec/len1 < args->min_overlap && isec/len2 < args->min_overlap ) continue;
                }

                // Process the comma-separated list of SRC annotations
                int nsrc_annots = 0;
                char **src_annots = hts_readlist(regitr_payload(args->itr,char*), 0, &nsrc_annots);

                int match = 1;
                if ( args->require_match )
                {
                    match = 0;
                    for (i=0; i<nsrc_annots; i++)
                    {
                        if ( !strcmp(dst_payload,src_annots[i]) ) { match = 1; break; }
                    }
                }
                
                if ( match )
                {
                    for (i=0; i<nsrc_annots; i++)
                    {
                        char *annot = args->default_hit ? args->default_hit : src_annots[i];
                        if ( args->only_unique && khash_str2int_has_key(dst_annots_hash,annot) ) continue;
                        if ( !first ) kputc(',',&str); else first = 0;
                        kputs(annot, &str);
                        if ( dst_annots_hash ) khash_str2int_set(dst_annots_hash, strdup(annot), 1);
                    }
                }

                for (i=0; i<nsrc_annots; i++) free(src_annots[i]);
                free(src_annots);
            }
        }
        chr_end[1] = tmp;
        if ( str.l )
        {
            if ( dst_payload && !args->new_column ) kputc(',',&line);
            else kputc('\t',&line);
            kputs(str.s, &line);
        }
        else if ( (!dst_payload && args->default_miss) || (args->new_column && args->default_miss) )
        {
            kputc('\t',&line);
            kputs(args->default_miss, &line);
        }
        printf("%s\n", line.s);

        // Clean up
        if ( dst_annots_hash )
            khash_str2int_destroy_free(dst_annots_hash);
        free(dst_payload);
    }

    hts_close(fp);
    free(str.s);
    free(line.s);
    regitr_destroy(args->itr);
    regidx_destroy(args->idx);
    free(args);

    return 0;
}

