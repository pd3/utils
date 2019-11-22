# Utilities
A random collection of tools, sorted in descendent order by the degree of usefulness they've had for me.

## [annot-regs](annot-regs.md)
An incredibly useful small program which allows fast and effortless execution of these common tasks:
 - find overlaps in two sets of genomic regions (e.g. two CNV callsets)
 - ... while matching records by a common column (e.g. sample name)
 - ... and transferring one or more fields from one file to another
 
See the [slides](https://raw.githubusercontent.com/pd3/utils/master/annot-regs/doc.annot-regs.pdf) for documentation and usage examples.

## [perm-test](perm-test.md)
An efficient permutation test for testing enrichment of calls in target vs background regions.
It is suitable especially for exomes because it avoids placing calls in inaccessible regions
while still correctly handling calls and regions of different sizes.

See the [slides](https://raw.githubusercontent.com/pd3/utils/master/perm-test/doc.perm-test.pdf) for documentation and usage examples.

## [dist](dist.md)
A memory-efficient program for collecting empirical distributions with very long tails.


