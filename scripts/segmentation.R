#!/usr/bin/env Rscript

## format
### file_non0
file_non0 <- function(f) { 
    sum(file.size(f), na.rm=T) > 0
}

### Unit Converter
unit_converter <- function(tx) { 
    div <- findInterval(as.numeric(gsub("\\,", "", tx)), 
                        c(0, 1e3, 1e6, 1e9, 1e12) )
    paste(round( as.numeric(gsub("\\,","",tx))/10^(3*(div-1)), 2), 
      c("","K","M","B","T")[div] )}

library(DNAcopy)
### Allele Frequency caller
segmentation = function(cr_dt, sampleid) {
    # remove 'chr'
    cr_dt[, CONTIG := gsub('^chr', '', CONTIG)]
    cr_dt[, COPY_RATIO := 2^LOG2_COPY_RATIO]

    cr_cna = CNA(genomdat = cr_dt$COPY_RATIO, 
            chrom = cr_dt$CONTIG,
            maploc = cr_dt$START,
            sampleid=sampleid, presorted=TRUE)
    cr_cna = smooth.CNA(cr_cna)
    fd = segment(cr_cna)

    fwrite(fd$output, 'EAC-M001-IMEAC1.seg', sep="\t")



    cn_dt[, minor_AF := min(H1_MEAN_AF, 1-H1_MEAN_AF), by = .(CONTIG, START, END)]

    cn_cna = CNA(genomdat = cn_dt$minor_AF, 
            chrom = cn_dt$CONTIG,
            maploc = cn_dt$START,
            sampleid=sampleid, presorted=TRUE)
    cn_cna = smooth.CNA(cn_cna)
    fd = segment(cn_cna)

    fwrite(fd$output, 'EAC-M001-IMEAC1.minor_AF.seg', sep="\t")



    }
    # generate a new column
    hets_dt[, POSITION2 := POSITION]
    setkey(hets_dt, CONTIG, POSITION, POSITION2)
    return(hets_dt)
}



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
    make_option(c("-o", "--output_prefix"), type = "character", default = "plotACN",
     help = "Output prefix [%default]", metavar = "character"))

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_ratio.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The first four columns of <copy_ratio.tsv> must be CONTIG, START, END and LOG2_COPY_RATIO.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

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
cr_dt = fread(r, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)

cn_dt = fread(n, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
cn_dt[, :fd]



## format
### process the columns in Copy Ratio table
r_dt[, CONTIG := gsub('^chr', '', CONTIG)]
r_dt = r_dt[, .(CONTIG, START, END, LOG2_COPY_RATIO)]
setkey(r_dt, CONTIG, START, END)







### process the columns in cytoBand table
b_dt[, chr := gsub('^chr', '', chr)]
setkey(b_dt, chr, start, end)
### call Allele Frequency
if (!file.exists(p)) p_dt = data.table()
if (!file.exists(e) | x == 1) e_dt = data.table()
if (nrow(e_dt) > 0) e_dt = AF_caller(hets_dt=e_dt, phased_dt=p_dt)
### infer threshold for hets number
if ((nrow(e_dt) > 0) & thn == -1) thn = floor((r_dt[1, END] - r_dt[1, START] + 1) / 25000 * 5)
### infer Allelic Depth
if ((nrow(e_dt) > 0) & ta == -1) ad = AD_caller(cr_dt=r_dt, hets_dt=e_dt, purity=u, ploidy=l, min_hets_num=thn)
if ((nrow(e_dt) == 0) & ta == -1) ad = 3
if (ta == -1) ta = ifelse(floor(ad / 2) < 3, 3, floor(ad / 2))
### clip Genomic Region
if (grepl('[0-9XY]+', g)) {
    g_dt = GR_extractor(g_region=g)
    r_dt = GR_clipper(data_dt=r_dt, ref_dt=g_dt)
} else {
    g = 'all'
    g_dt = data.table()
}
### clip Gene Table
if (exists('gene_dt')) {
    gene_dt = GT_extractor(data_dt=r_dt, gene_dt=gene_dt)
} else {
    gene_dt = data.table()
}


## calculate
### call Copy Number
cn_dt = CN_caller(cr_dt=r_dt, hets_dt=e_dt, gene_dt=gene_dt, purity=u, ploidy=l, min_hets_num=thn, neu_cn=x)
### filter Structural Variants
if (file.exists(n)) sv_dt = SV_filter(sv_dt=sv_dt, sample_id=s, min_sv_AD=ta, min_sv_span=tp)


## map
### assemble contig
if (!file.exists(d)) contig_dict_df = data.table()
contig_dt = contig_assembler(ref_dict_df=contig_dict_df, g_region_dt=g_dt, cn_dt=cn_dt)
### map cn_dt to contig
cn_gdt = cn2contig(cn_dt=cn_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### map sv_dt to contig
if (!file.exists(n)) sv_gdt = data.table()
if (file.exists(n)) sv_gdt = sv2contig(sv_dt=sv_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### map cytoband to contig
if (!file.exists(b)) b_gdt = data.table()
if (file.exists(b)) b_gdt = cb2contig(b_dt=b_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### convert to data.frame
cn_df = CN_format(cn_dt, e_dt, p)
g_df = as.data.frame(g_dt)
if (file.exists(n)) sv_df = as.data.frame(sv_dt)
contig_df = as.data.frame(contig_dt)
cn_gdf = as.data.frame(cn_gdt)
sv_gdf = as.data.frame(sv_gdt)
b_gdf = as.data.frame(b_gdt)


## plot
### print scatter plot for CN
pdf(paste0(o, '.CN_scatter.pdf'), pw, ph, useDingbats=FALSE)
    cs_gb=CN_scatter(cn_gdf=cn_gdf, sv_gdf=sv_gdf, contig_df=contig_df)
    grid.draw(cs_gb)
dev.off()
### print scatter plot for AF
if (nrow(e_dt) > 0) {
    if(file.exists(p)) {
        af_dt = copy(cn_dt)
        af_dt[, `:=`(XAF1 = MEAN_XAF, XAF2 = 1 - MEAN_XAF), ]
    } else {
        af_dt = copy(e_dt)
        af_dt[, `:=`(XAF1 = XAF, XAF2 = 1 - XAF, MIDDLE = POSITION), ]
    }
    af_dt = af_dt[, .(CONTIG, MIDDLE, XAF1, XAF2)]
    af_dt = sAF_caller(af_dt)
    af_gdt = af2contig(cn_dt=af_dt, contig_dt=contig_dt, g_region_dt=g_dt)
    pdf(paste0(o, '.AF_scatter.pdf'), pw, ph, useDingbats=FALSE)
        as_gb = AF_scatter(af_gdt=af_gdt, sv_gdf=sv_gdf, contig_df=contig_df)
        grid.draw(as_gb)
    dev.off()
}
### print contour plot
if (nrow(e_dt) > 0) {
    cn4_max = 2^ceiling(log2(max(cn_df[, 4])))
    cn5_max = 2^ceiling(log2(max(cn_df[, 5])))
    pdf(paste0(o, '.CN_contour.pdf'), useDingbats=FALSE, 
      width = (3 + (ceiling(log2(cn4_max))-1)/2*2) / 2.54 * 1.5, 
      height = (3 + (ceiling(log2(cn5_max))-1)/2*2) / 2.54 * 1.5)
        cc_gb=CN_scatter2(cn_df=cn_df, xmax=cn4_max, ymax=cn5_max)
        grid.draw(cc_gb)
    dev.off()
}


## write
### print tables
write.table(cn_df, paste0(o, '.CN.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
if (nrow(sv_gdf) > 0) write.table(sv_gdf, paste0(o, '.SV.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
cat(paste('sample', 'region', 'min_allelic_depth', 'min_hets_num', sep='\t'), '\n', sep='')
cat(paste(s, g, ta, thn, sep='\t'), '\n', sep='')



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
