### Calculate observed 1D-folded-SFS (minor allele frequency spectrum) for stairway plot
### Use Stacks output (populations.log and vcf files).
### Missing data was compensated by a bootstrapping.
### Output basic information and SFS as .csv and .txt files, respectively.

calc.1d.msfs <- function (infile.dir, outfile.basic.info, outifle.sfs) {
### Extract information from log file of populations
    logfile <- paste(infile.dir, "populations.log", sep = "")
    d.log <- readLines(con = logfile, warn = FALSE)
    temp <- unlist(strsplit(grep("population map contained", d.log, value = TRUE),
                            "population map contained "))
    n.ind <- as.numeric(strsplit(temp[2], " ")[[1]][1]) # Numner of individuals
    
    temp <- unlist(strsplit(grep("Kept", d.log, value = TRUE), " "))
    n.stacks.locus <- as.numeric(temp[grep("Kept", temp) + 1]) # Number of stacks-loci including 0 SNP loci
    n.all.site <- as.numeric(temp[grep("composed", temp) + 2]) # Number of all sites including monomorphic sites
    n.snp <- as.numeric(temp[grep("filtered,", temp) + 1]) # Number of SNPs (variant sites)
    mean.locus.length <- n.all.site/n.stacks.locus
    mean.snp.per.locus <- n.snp/n.stacks.locus
    
    ## output
    out <- data.frame(Variable = c("n.ind", "n.stacks.locus", "n.all.site",
                                   "n.snp", "mean.locus.length",
                                   "mean.snp.per.locus"),
                      Value = c(n.ind, n.stacks.locus, n.all.site, n.snp,
                                mean.locus.length, mean.snp.per.locus))
    write.table(out, file = outfile.basic.info, sep = ",", quote = FALSE,
                row.names = FALSE)
    rm(out)

### Print locus information
    print(paste("Number of individuals is", n.ind))
    print(paste("Total number of stacks-loci is", n.stacks.locus))
    print(paste("Total number of sites is", n.all.site))
    print(paste("Total number of SNP is", n.snp))
    print(paste("Average length of stack-loci is", mean.locus.length))
    print(paste("Average number of SNP per stacks-loci is", mean.snp.per.locus))

### Read VCF, compensate missing and calculate SFS
    vcffile <- paste(infile.dir, "populations.snps.vcf", sep = "")
    d <- read.table(file = vcffile, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
    names(d) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT"),
                  paste("Ind", 1:n.ind, sep = ""))

    out <- data.frame(NMA = 0:n.ind, # Number of minor allele
                      Count = rep(0, n.ind+1))
    
    out$Count[1] <- n.all.site - n.snp # Number of monomorphic sites
    
    for(i in 1:nrow(d)) {
        temp <- d[i, grep("Ind", names(d), value = TRUE)]
        temp <- substring(temp, 1, 3)
        allele <- unlist(strsplit(temp, "/"))
        if(any(allele == ".")) { # When there are any missing data, complement by bootstrapping
            n.missing <- sum(allele == ".") # Number of missing allele
            allele <- allele[allele != "."] # Remove missing
            allele <- c(allele,
                        sample(allele, n.missing, replace = TRUE))
        }
        NMA <- sort(table(allele))[1] # Number of minor allele
        out$Count[out$NMA == NMA] <- out$Count[out$NMA == NMA] + 1
    }
    
    header <- paste("d0_", out$NMA, sep = "")

    sink(outfile.sfs)
    cat(header, sep = "\t")
    cat("\n")
    cat(out$Count, sep = "\t")
    cat("\n")
    sink()
}

#set.seed(46) # Fix random seed
#infile.dir <- "~/works/works16/ibuki/data/stacks/210813VS_Hokkaido/analysis01/"
#outfile.basic.info <- "basic_info_VS10inds.csv"
#outfile.sfs <- "observed_SFS_VS10inds.txt"
#calc.1d.msfs (infile.dir, outfile.basic.info, outfile.sfs)

