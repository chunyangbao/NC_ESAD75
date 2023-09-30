#!/usr/bin/env Rscript

source('~/bin/ACN_fun.R')

# libraries
## basic
library(optparse)
library(data.table)
## ggplot
library(ggplot2)
library(ggforce)
library(grid)
library(gridExtra)
library(gtable)
library(ggrepel)


# options
V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table, ggplot2, ggforce, grid, gridExtra, gtable"

a = commandArgs(trailingOnly = TRUE)

option_list = list(
    make_option(c("-u", "--purity"), type = "double", default = "1",
     help = "Purity [%default]", metavar = "double"),
    make_option(c("-l", "--ploidy"), type = "double", default = "2",
     help = "Ploidy [%default]", metavar = "double"),
    make_option(c("--plot_width"), type = "double", default = "4",
     help = "Width of the ouput plot [%default]", metavar = "double"),
    make_option(c("--plot_height"), type = "double", default = "4",
     help = "Height of the ouput plot [%default]", metavar = "double"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "plotACN",
     help = "Output prefix [%default]", metavar = "character"))

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_ratio.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The first five columns of <copy_ratio.tsv> must be CONTIG, START, END, MAJOR_CN and MINOR_CN.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

u = opt$options$purity
l = opt$options$ploidy
pw = opt$options$plot_width
ph = opt$options$plot_height
o = opt$options$output_prefix

r = opt$args[1]


# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    r = ifelse(r == '-', 'file:///dev/stdin', r)
} else {
    r = ifelse(r == '-', 'cat /dev/stdin', r)
}
### fread
cn_dt = fread(r, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
cn_df = as.data.frame(cn_dt)

## plot
### print contour plot
pdf(paste0(o, '.CN_scatter.pdf'), pw, ph, useDingbats=FALSE)
    cc_gb=CN_scatter2(cn_df=cn_df, purity=u, ploidy=l)
    grid.draw(cc_gb)
dev.off()
