Calculate 1D-folded SFS from VCF and population.log files of stacks using R
======
* This script calculate observed 1D-folded SFS (minor allele frequency spectrum) for stairway plot.
* It uses stacks output (populations.log and vcf files).
* Missing data is compensated by a bootstrapping.
* It outputs basic information and SFS as .csv and .txt files, respectively.


Example
------
Example for calculation of 1D folded SFS

```R
    source("r230327calc_sfs.r")
  
    set.seed(46) # Fix random seed
    infile.dir <- "~/works/works16/ibuki/data/stacks/210813VS_Hokkaido/analysis01/" # Path for stacks output directory
    outfile.basic.info <- "basic_info_VS10inds.csv" # output basic information file name
    outfile.sfs <- "observed_SFS_VS10inds.txt" # output 
  
    calc.1d.msfs (infile.dir, outfile.basic.info, outfile.sfs)
```
