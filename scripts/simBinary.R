#!/usr/bin/env Rscript

# simBinary
simBinary <- function(bin_ma) {
    col_comb <- combn(colnames(bin_ma), 2)
    sim_df <- data.frame(sample_1=col_comb[1, 1], 
                         sample_2=col_comb[1, 1],
                         similarity=1)
    for (i in 1:ncol(col_comb)) {
        if(sum(grepl(col_comb[1, i], sim_df[, 1]))==0) {
            i0_df <- data.frame(sample_1=col_comb[1, i],
                                sample_2=col_comb[1, i],
                                similarity=1)
            sim_df <- rbind(sim_df, i0_df)
        }
        i_ma <- bin_ma[,c(col_comb[1, i], col_comb[2, i])]
        nHits_both <- sum(rowSums(i_ma) == 2)
        nHits_1 <- sum(i_ma[, 1])
        nHits_2 <- sum(i_ma[, 2])
        nHits <- nHits_1 + nHits_2 - nHits_both
        i1_df <- data.frame(sample_1=col_comb[1, i],
                            sample_2=col_comb[2, i],
                            similarity=as.numeric(nHits_both / nHits))
        i2_df <- data.frame(sample_1=col_comb[2, i],
                            sample_2=col_comb[1, i],
                            similarity=as.numeric(nHits_both / nHits))
        i_df <- rbind(i1_df, i2_df)
        sim_df <- rbind(sim_df, i_df)
    }
    return(sim_df)
}

 


library("optparse")
library("data.table")

library("ComplexHeatmap")

V <- "Version: 1.0"
D <- "Depends: R (>= 3.1.0), optparse, data.table, ComplexHeatmap"

a <- commandArgs(trailingOnly=TRUE)

option_list <- list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv>\n\tBy default <input.tsv> is '-' (stdin).", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

# options
o <- opt$options$output


# input
## R version
rv0 <- as.numeric(R.version$major)
rv1 <- as.numeric(R.version$minor)
## stdin
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1])
} else {
    d <- ifelse(opt$args[1]=='-', 'cat /dev/stdin', opt$args[1])
}

# read
d <- fread(d, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
d_ma <- as.matrix(d[,-1], rownames = as.character(d[[1]]))

# calculate similarity
sim_df <- simBinary(d_ma)

# df to ma
sim_df <- setorder(sim_df, sample_1, sample_2)
sim_df <- dcast(sim_df, sample_1 ~ sample_2, value.var = 'similarity')
sim_ma <- as.matrix(sim_df[,-1])
rownames(sim_ma) <- as.character(sim_df[[1]])
sim_ma[is.na(sim_ma)] <- 1

# plot
col_fun = circlize::colorRamp2(c(0, 1), c("white", "black"))
pdf(paste0(o, ".similarity.pdf"), width = 11, height = 11, , useDingbats=FALSE)
    Heatmap(sim_ma, name = "similarity", col = col_fun, rect_gp = gpar(type = "none"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = NA, fill = NA))
            grid.circle(x = x, y = y, r = width * 0.3, 
                gp = gpar(fill = col_fun(sim_ma[i, j]), col = NA))
        }, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE)
dev.off()

fwrite(sim_df, paste0(o, ".similarity.tsv"), sep = "\t")
