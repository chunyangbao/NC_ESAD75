#!/usr/bin/env Rscript

# functions
## filter
### Flip by Reference per Bin
FlipbyRefperBin = function(gt) {
    fgt_dt = copy(gt) # flipped GT
    fgt_dt[, GT_diff := GT_e * GT_r]
    fgt_dt[, GT_sum := sum(GT_diff, na.rm = TRUE), by = .(CHROM, START_g, END_g)]
    fgt_dt[, GT_N := sum(abs(GT_diff), na.rm = TRUE), by = .(CHROM, START_g, END_g)]
    fgt_dt[, GT_ratio := GT_sum/GT_N]
    fgt_dt[, GT_ratio := round(GT_ratio)]
    fgt_dt[, GT_ratio := ifelse(GT_ratio == -1, -1, 1)]
    fgt_dt[, GT_ratio := ifelse(is.na(GT_ratio), 1, GT_ratio)]
    fgt_dt[, GT_e := GT_e * GT_ratio]
    fgt_dt = fgt_dt[, .(CHROM, START, END, REF, ALT, GT_e)]
    return(fgt_dt)
}
### Flip by Reference
FlipbyRef = function(gt) {
    fgt_dt = copy(gt) # flipped GT
    fgt_dt[, GT_diff := GT_e * GT_r]
    # iterate each chr
    for (i in unique(gt[,CHROM])) {
        repeat {
            gt_i = fgt_dt[CHROM == i]
            ns_i = nrow(gt_i[GT_diff == -1]) # number of switching errors
            if (ns_i == 0) break
            ps1_i = gt_i[GT_diff == -1, START][1] # position of the 1st switching error
            fgt_dt[(CHROM == i) & (START >= ps1_i), GT_e := GT_e * -1]
        }
    }
    fgt_dt = fgt_dt[, .(CHROM, START, END, REF, ALT, GT_e)]
    return(fgt_dt)
}


library(optparse)
library(data.table)

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table, vcfR, e1071"

a = commandArgs(trailingOnly = TRUE)

option_list = list(
    make_option(c("-g", "--ref_region"), type="character", default="",
     help="Reference SNPs in VCF format [%default].", metavar="character"),
    make_option(c("-n", "--ref_intervals"), type="character", default="",
     help="Reference intervals in GATK4 interval_list format [%default].", metavar="character"),
    make_option(c("-r", "--ref_vcf"), type="character", default="",
     help="Reference SNPs in VCF format [%default].", metavar="character")
)

opt_parser = OptionParser(usage = "usage: %prog [options] <input.vcf>\n\t
    This program flips fliping vector accroding to reference GT in a VCF \n\t
    and output a flipped VCF. This required input file should be in VCF \n\t
    format (separator: '\\t', stdin: '-')", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

g = opt$options$ref_region
n = opt$options$ref_intervals
r = opt$options$ref_vcf

e = opt$args[1]

# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    e = ifelse(e == '-', 'file:///dev/stdin', e)
    if(grepl('\\.gz$', e)) e_dt = fread(cmd = paste0('zcat ', e), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
    if(file.exists(r)) {if(grepl('\\.gz$', r)) r_dt = fread(cmd = paste0('zcat ', r), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
} else {
    e = ifelse(e == '-', 'cat /dev/stdin', e)
    if(grepl('\\.gz$', e)) e_dt = fread(paste0('zcat ', e), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
    if (file.exists(r)) {if(grepl('\\.gz$', r)) r_dt = fread(paste0('zcat ', r), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
}
### fread
if(file.exists(n)) {n_dt = fread(n, sep = "\t", select = c(1,2,3), fill = TRUE, colClasses=list(character=1), stringsAsFactors = FALSE, header = FALSE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
if(!grepl('\\.gz$', e)) e_dt = fread(e, skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(r)) {if(!grepl('\\.gz$', r)) r_dt = fread(r, skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}

## format
if (!file.exists(n)) n_dt = data.table()
if (nrow(n_dt) > 0) {
    colnames(n_dt) = c('CHROM', 'START_g', 'END_g')
    n_dt = n_dt[!CHROM %like% "@"]
    n_dt[, START_g := as.numeric(START_g)]
    n_dt[, END_g := as.numeric(END_g)]
    setkey(n_dt, CHROM, START_g, END_g)
}
### update column names for 'CHROM' column
colnames(e_dt) = sub('#', '', colnames(e_dt))
colnames(r_dt) = sub('#', '', colnames(r_dt))
### update column names for 'GT' column
colnames(e_dt)[grep("|", e_dt[1], fixed=T)] = "GT_e"
colnames(r_dt)[grep("|", r_dt[1], fixed=T)] = "GT_r"
### subset of hets sites
e_dt = e_dt[(GT_e == "0|1")|(GT_e == "1|0")]
r_dt = r_dt[(GT_r == "0|1")|(GT_r == "1|0")]
### update colnames
e_dt[, START := POS]; e_dt[, END := POS]
r_dt[, START := POS]; r_dt[, END := POS]
### select columns
e_dt = e_dt[,.(CHROM, START, END, REF, ALT, GT_e)]
r_dt = r_dt[,.(CHROM, START, END, REF, ALT, GT_r)]
### select columns
setkey(e_dt, CHROM, START, END, REF, ALT)
setkey(r_dt, CHROM, START, END, REF, ALT)

## process
### merge
re_dt = merge(e_dt, r_dt, all.x=TRUE)
re_dt = re_dt[(!is.na(GT_e))&(!is.na(GT_r))]
setorder(re_dt, CHROM, START)
re_dt[, GT_e := gsub("0|1", -1, GT_e, fixed=T)]
re_dt[, GT_e := gsub("1|0", 1, GT_e, fixed=T)]
re_dt[, GT_r := gsub("0|1", -1, GT_r, fixed=T)]
re_dt[, GT_r := gsub("1|0", 1, GT_r, fixed=T)]
re_dt[, GT_e := as.numeric(GT_e)]
re_dt[, GT_r := as.numeric(GT_r)]
setkey(re_dt, CHROM, START, END)
### flip
if (nrow(n_dt) > 0) {
    re_dt = foverlaps(re_dt, n_dt)
    re1_dt = FlipbyRefperBin(re_dt)
} else {
    re1_dt = FlipbyRef(re_dt)
}

### format
re1_dt[, POS := START]
re1_dt = re1_dt[, .(CHROM, POS, REF, ALT, GT_e)]
re1_dt[, GT_e := as.character(GT_e)]
re1_dt[GT_e == "1", GT_e := "1|0"]
re1_dt[GT_e == "-1", GT_e := "0|1"]
colnames(re1_dt)[1] = "#CHROM"
## write
fwrite(re1_dt, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)

### THE END ###
