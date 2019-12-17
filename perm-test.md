```
License: The MIT/Expat license
About: Run enrichment (permutation) test. It is efficient, does not attempt to place calls
       in inaccessible regions, and correctly handles calls and regions of different sizes.
       All coordinates are 1-based, inclusive, unless the file names have ".bed", ".bed.gz"
       or ".bed.bgz" extension.
Usage: perm-test [OPTIONS]
Options:
   -b, --background-regs FILE      Background regions, expected not to be enriched: chr,beg,end
   -c, --calls FILE                Calls: chr,beg,end
   -d, --debug-regions             Print the spliced regions (and stop)
   -f, --ref-fai FILE              Chromosome lengths, given for example as fai index: chr,length
       --no-bg-overlap             Permuted variants must not overlap background regions
   -o, --output FILE               Place output in FILE
   -m, --max-call-length INT       Skip big calls [10e6]
   -n, --niter NUM1[,NUM2]         Number of iterations: total (NUM1) and per-batch (NUM2, reduces memory)
   -s, --random-seed INT           Random seed
   -t, --target-regs FILE          Target regions, expected to be enriched: chr,beg,end

Examples:
   # Take a set of "background" regions (e.g. exome baits), "target" regions (e.g. a set of genes)
   # and calculate how likely it is for a set of calls (e.g. CNVs) to overlap the target regions given
   # the accessible background+target regions. The program ensures that overlapping and duplicate
   # regions are handled correctly.
   perm-test -b baits.txt -t genes.txt -c calls.txt -f ref.fai -n 1e9,1e8

   # Same as above but this time count only permuted variants that overlap target regions but do not
   # overlap any background region (e.g. CNVs that hit only a non-coding region but not a coding sequence).
   # Note that target regions take precedence over background regions in case the sets overlap. If that is
   # not desired, use `bedtools subtract -a noncoding.bed -b coding.bed` to clip and write regions from
   # file A not present in B.
   perm-test -b coding.bed -t noncoding.bed -c calls.txt -f ref.fai -n 1e9,1e8 --no-bg-overlap

```
