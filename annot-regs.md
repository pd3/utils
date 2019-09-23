```
About: Annotate regions in DST file with texts from overlapping regions in SRC file.
       The transfer of annotations can be conditioned on matching values in one or more
       columns (-m), multiple columns can be transferred (-t).
       All indexes and coordinates are 1-based and inclusive.
Usage: annot-regs [OPTIONS] DST
Options:
       --allow-dups                Add annotations multiple times
   -a, --annotate list             Add special annotations:
                                       frac .. fraction of the destination region with an overlap
                                       nbp  .. number of source base pairs in the overlap
   -c, --core src:dst              Core columns [chr,beg,end:chr,beg,end]
   -d, --dst-file file             Destination file
   -H, --ignore-headers            Use numeric indexes, ignore the headers completely
   -m, --match src:dst             Require match in these columns
   -o, --overlap float             Minimum required overlap (non-reciprocal, unless -r is given)
   -r, --reciprocal                Require reciprocal overlap
   -s, --src-file file             Source file
   -t, --transfer src:dst          Columns to transfer. If src column does not exist, interpret
                                   as the default value to use. If the dst column does not exist,
                                   a new column is created.
       --version                   Print version string and exit
Examples:
   # Header is present, match and transfer by column name
   annot-regs -s src.txt.gz -d dst.txt.gz -c chr,beg,end:chr,beg,end -m type,sample:type,smpl -t tp/fp:tp/fp

   # Header is not present, match and transfer by column index (1-based)
   annot-regs -s src.txt.gz -d dst.txt.gz -c 1,2,3:1,2,3 -m 4,5:4,5 -t 6:6

   # If the dst part is not given, the program assumes that the src:dst columns are identical
   annot-regs -s src.txt.gz -d dst.txt.gz -c chr,beg,end -m type,sample -t tp/fp
```
