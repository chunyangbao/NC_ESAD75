#!/usr/bin/env Rscript

source('~/bin/ACN_fun.R')

# libraries
## basic
library(optparse)
library(data.table)


# options
V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table, ggplot2, ggforce, grid, gridExtra, gtable"

a = commandArgs(trailingOnly = TRUE)

option_list = list(
    make_option(c("-e", "--hets_tsv"), type = "character", default = "",
     help = "Het sites data file in TSV format. The first four columns must be \n\t\tCONTIG, POSITION, REF_COUNT and ALT_COUNT. (optinal) [%default]", metavar = "character"),
    make_option(c("-u", "--purity"), type = "double", default = "1",
     help = "Purity [%default]", metavar = "double"),
    make_option(c("-l", "--ploidy"), type = "double", default = "2",
     help = "Ploidy [%default]", metavar = "double"),
    make_option(c("-x", "--neu_cn"), type = "integer", default = "2",
     help = "Copy number of X chromosome [%default] \n\t\te.g. 2 for female's X, 1 for male's X", metavar = "integer"),
    make_option(c("--th_hetn"), type = "double", default = "-1",
     help = "Threshold for hets number (roughly, 0 for 5kb, 2 for 10kb, 5 for 25kb) [%default]", metavar = "double")
)

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_ratio.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The first four columns of <copy_ratio.tsv> must be CONTIG, START, END and LOG2_COPY_RATIO.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

e = opt$options$hets_tsv
u = opt$options$purity
l = opt$options$ploidy
x = opt$options$neu_cn
thn = opt$options$th_hetn

r = opt$args[1]


# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    r = ifelse(r == '-', 'file:///dev/stdin', r)
    e = ifelse(e == '-', 'file:///dev/stdin', e)
} else {
    r = ifelse(r == '-', 'cat /dev/stdin', r)
    e = ifelse(e == '-', 'cat /dev/stdin', e)
}
### fread
r_dt = fread(r, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(e)) e_dt = fread(e, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)


## format
### process the columns in Copy Ratio table
r_dt[, CONTIG := gsub('^chr', '', CONTIG)]
r_dt = r_dt[, .(CONTIG, START, END, LOG2_COPY_RATIO)]
setkey(r_dt, CONTIG, START, END)

## calculate
### call Copy Number
cn_dt = CN_caller(cr_dt=r_dt, hets_dt=e_dt, purity=u, ploidy=l, min_hets_num=thn, neu_cn=x)

## write
### print tables
write.table(cn_dt, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)



# test
#cn_dt = data.table(CONTIG=c(rep(1,5),2,2,2,3,3),
#    START=c(1,2,6,10,12,18,3,11,2,7),END=c(3,5,8,11,14,22,8,17,6,12),
#    MAJOR_CN=c(0,0.8,1,1.1,1.8,2,2.2,4,22,54), MINOR_CN=c(0,0.4,0.8,1,0.86,0.9,0.85,0.76,11,20))
#sv_dt = data.table(CONTIG=c(rep(1,4),2,2,2,2,3,3),
#    POS=c(2,3,5,10,5,7,8,11,5,7),POS2=c(2,6,5,11,12,9,15,17,12,9),
#    SV_type=c('DEL','AMP','AMP','DEL','DEL','INV','DEL','AMP','AMP','INV'))
#d=''

#s = 'SM_9_11_EAC1'
#r = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p1p2/denoised_copy_ratios_tumor/9_11_EAC1.denoisedCR.tsv'
#e = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p1p2/allelic_counts_tumor/9_11_EAC1.allelicCounts.tsv'
#n = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/SV_evolution/CZ_bkps/EAC_9_11_bkps.txt'
#d = '/cga/bass/Chunyang/ref/hg19/Homo_sapiens_assembly19.1_22.dict'
#k = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/PhylogicNDT/Mutect2_gga_filtered1GF_forecalled_census_bed/EAC-9_11-EAC1.census.mut_cnv.bed'
#p = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p1p2/het_allelic_counts_normal_t0_phased1f_vcf/EAC-9_11-NORM.hets.phased1f.vcf.gz'
#b = '/cga/bass/Chunyang/ref/hg19/cytoBand.txt'
#u = 0.4
#l = 3.85
#o = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p1p2/allelic_counts_tumor_plotACN/EAC-9_11-EAC1.chr17.5kb'
#g = 'chr17:35000000-45000000'
#thn = 0

### scale for x axis (developing)
#x_max = contig_df[nrow(contig_df), 'ends']
#x_max = x_max + 2 * x_max * .01
#if ((nrow(contig_df) == 1)&(ncol(g_df) == 1)) {
#    gs_width = x_max / 40000000
#} else {
#    gs_width = 7
#}
### scale for y axis
#if (nrow(sv_gdf) > 0) {
#    gs_height = (6 + (3 + (ceiling(log2(cn_max))-1)/2*2)) / 2 / 2.54
#} else {
#    gs_height = (0 + (3 + (ceiling(log2(cn_max))-1)/2*2)) / 2 / 2.54
#}
