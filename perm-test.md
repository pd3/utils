```
About: Run enrichment (permutation) test. It is efficient, does not attempt to place calls
       in inaccessible regions, and correctly handles calls and regions of different sizes.
Usage: perm-test [OPTIONS]
Options:
   -b, --background-regs FILE      Background regions, expected not to be enriched: chr,beg,end
   -c, --calls FILE                Calls: chr,beg,end
   -d, --debug-regions             Print the spliced regions (and stop)
   -f, --ref-fai FILE              Chromosome lengths, given for example as fai index: chr,length
   -o, --output FILE               Place output in FILE
   -m, --max-call-length INT       Skip big calls [10e6]
   -n, --niter NUM1[,NUM2]         Number of iterations: total (NUM1) and per-batch (NUM2, reduces memory)
   -s, --random-seed INT           Random seed
   -t, --target-regs FILE          Target regions, expected to be enriched: chr,beg,end
Example:
   perm-test -b bg.txt -t tgt.txt -c calls.txt -f ref.fai -n 1e9,1e8
```
