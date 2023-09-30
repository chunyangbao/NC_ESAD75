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
    make_option(c("-e", "--hets_tsv"), type = "character", default = "",
     help = "Phased genotypes in VCF format (single sample VCF). Please refer to 'https://imputation.sanger.ac.uk/' (optinal) [%default]", metavar = "character"),
    make_option(c("-u", "--purity"), type = "double", default = "1",
     help = "Purity [%default]", metavar = "double"),
    make_option(c("-l", "--ploidy"), type = "double", default = "2",
     help = "Ploidy [%default]", metavar = "double"),
    make_option(c("--th_hetn"), type = "double", default = "-1",
     help = "Threshold for hets number (roughly, 0 for 5kb, 2 for 10kb, 5 for 25kb) [%default]", metavar = "double"),
    make_option(c("-x", "--neu_cn"), type = "integer", default = "2",
     help = "Copy number of X chromosome [%default] \n\t\te.g. 2 for female's X, 1 for male's X", metavar = "integer")
    )

opt_parser <- OptionParser(usage = "usage: %prog [options] <allelic_counts.tsv>\n
    the input files should be in TSV format (separator: '\\t', stdin: '-')\n
    <genotype.vcf> is a VCF-like file and must include the following columns: CHROM, POS, REF, ALT, GT.\n 
    <allelic_counts.tsv> is a Allelic-counts file generated from GATK4 CollectAllelicCounts.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt <- parse_args(opt_parser, args = a, positional_arguments = TRUE)

e <- opt$options$hets_tsv
u <- opt$options$purity
l <- opt$options$ploidy
x = opt$options$neu_cn
thn = opt$options$th_hetn

d <- opt$args[1]

# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    d = ifelse(d == '-', 'file:///dev/stdin', d)
} else if (rv0>=4){
    d = ifelse(d == '-', 'file:///dev/stdin', d)
} else {
    d = ifelse(d == '-', 'cat /dev/stdin', d)
}
### fread
data = fread(d, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(e)) hets <- fread(e, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (!file.exists(e)) hets <- NULL

if ((nrow(hets) > 0) & thn == -1) thn = floor((data[1, END] - data[1, START] + 1) / 25000 * 5)

## format
### call Allele Frequency
cn = CN_caller(data, hets = hets, purity = u, ploidy = l, min_hets = thn, neu_cn = x)

## write
### print tables
fwrite(cn, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
