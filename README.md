# TLA
Telomere local assembly

# Install   
TLA relies on [htslib](https://github.com/samtools/htslib) and [spoa](https://github.com/rvaser/spoa)
```
make HTSLIB_PATH=/path/to/htslib SPOA_PATH=/path/to/spoa
```

# Usage
```
TLA: Telomere local assembly   
Usage: TLA [options]   
options:   
-b string     alignment file in bam format   
-m int        minimum length of contigs will be used (default:500kb)   
-d int        the distance to the start or end of contigs, used to extract reads in the BAM file (default:10kb)   
-s int        minimum length of perfect tandem repeats region(default:100bp)   
-p int        minimum percentage of reads contained perfect tandem repeats (default:0.3)   
-c int        coverage threshold for the generation of consensus sequences (default:2)   
```

![alt text](https://github.com/Wenfei-Xian/TLA/blob/main/TLA_github.png)
