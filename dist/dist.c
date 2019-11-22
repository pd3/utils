/* 
    Copyright (C) 2019 Genome Research Ltd.
    
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
#include <getopt.h>
#include <assert.h>
#include <stdarg.h>
#include "dist.h"
#include "version.h"

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Create an empirical distribution of positive integer values with a very long tail\n"
        "        using a logarithmic binning to keep memory low\n"
        "Usage: dist [OPTIONS] DST\n"
        "Options:\n"
        "   -n, --nprecise int              Number of orders of magnitude to represent exactly [4]\n"
        "   -v, --version                   Print version string and exit\n"
        "Examples:\n"
        "   # Represent values up to 1000 exactly, then use logarithmic binning\n"
        "   cat dat.txt | dist -n 3 \n"
        "\n"
        "   # Simple test to demonstrate the usage\n"
        "   for i in $(seq 1 50); do echo $i; done | dist -n 1\n"
        "\n";
}

int main(int argc, char **argv)
{
    int nprecise = 4;
    static struct option loptions[] =
    {
        {"version",no_argument,NULL,'v'},
        {"nprecise",required_argument,NULL,'n'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int c;
    while ((c = getopt_long(argc, argv, "hvn:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'v': printf("%s\n",UTILS_VERSION); return 0; break;
            case 'n': 
                nprecise = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse --nprecise %s\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%sVersion: %s\n\n", usage_text(), UTILS_VERSION); break;
        }
    }

    if ( isatty(fileno((FILE *)stdin)) ) error("%sVersion: %s\n\n", usage_text(), UTILS_VERSION);

    dist_t *dist = dist_init(nprecise);

    uint32_t num;
    while ( !feof(stdin) && fscanf(stdin,"%u",&num)==1 )
    {
        dist_insert(dist,num);
    }

    printf("#[beg\tend)\tcount\tdensity\n");

    uint32_t beg, end, i, n = dist_n(dist);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(dist, i, &beg, &end);
        if ( !cnt ) continue;

        // Print the interval, count and density
        printf("%u\t%u\t%"PRIu64"\t%f\n", beg, end, cnt, (double)cnt/(end-beg));
    }

    dist_destroy(dist);

    return 0;
}

