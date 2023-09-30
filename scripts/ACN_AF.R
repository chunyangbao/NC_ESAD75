#!/usr/bin/env Rscript

source('~/bin/ACN_fun.R')

# libraries
## basic
library(optparse)
library(data.table)


# options
V <- "Version: 1.0"
D <- "Depends: R (>= 3.4.0), optparse, data.table"

a <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("-g", "--genotype_vcf"), type = "character", default = "",
     help = "Phased genotypes in VCF format (single sample VCF). Please refer to 'https://imputation.sanger.ac.uk/' (optinal) [%default]", metavar = "character")
    )

opt_parser <- OptionParser(usage = "usage: %prog [options] <allelic_counts.tsv>\n
    the input files should be in TSV format (separator: '\\t', stdin: '-')\n
    <genotype.vcf> is a VCF-like file and must include the following columns: CHROM, POS, REF, ALT, GT.\n 
    <allelic_counts.tsv> is a Allelic-counts file generated from GATK4 CollectAllelicCounts.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt <- parse_args(opt_parser, args = a, positional_arguments = TRUE)

g <- opt$options$genotype_vcf

d <- opt$args[1]

# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0 < 3)|((rv0 == 3)&(rv1 < 5))) {
    d = ifelse(d == '-', 'file:///dev/stdin', d)
    if ((file.exists(g)) & (grepl('\\.gz$', g))) genotype = fread(cmd = paste0('zcat ', g), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
} else if (rv0 >= 4){
    d = ifelse(d == '-', 'file:///dev/stdin', d)
    if ((file.exists(g)) & (grepl('\\.gz$', g))) genotype = fread(cmd = paste0('zcat ', g), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
} else {
    d = ifelse(d == '-', 'cat /dev/stdin', d)
    if ((file.exists(g)) & (grepl('\\.gz$', g))) genotype = fread(paste0('zcat ', g), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
}
### fread
data <- fread(d, skip = 'CONTIG', sep = "\t", colClasses = list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if ((file.exists(g)) & (!grepl('\\.gz$', g))) genotype <- fread(g, skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (!file.exists(g)) genotype <- NULL

## format
### call Allele Frequency
af = AF_caller(data, genotype = genotype)

## write
### print tables
fwrite(af, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
