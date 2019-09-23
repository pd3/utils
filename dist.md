```
About: Create an empirical distribution of positive integer values with a very long tail
        using a logarithmic binning to keep memory low
Usage: dist [OPTIONS] DST
Options:
   -n, --nprecise int              Number of orders of magnitude to represent exactly [4]
   -v, --version                   Print version string and exit
Examples:
   # Represent values up to 1000 exactly, then use logarithmic binning
   cat dat.txt | dist -n 3 

   # Simple test to demonstrate the usage
   for i in $(seq 1 50); do echo $i; done | ./dist -n 1
```
