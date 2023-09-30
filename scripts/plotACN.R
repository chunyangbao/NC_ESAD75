#!/usr/bin/env Rscript

# functions
## color
### color vector
map = NULL
### Copy Number colors
map['MAJOR_CN'] = 'darkcyan'
map['MINOR_CN'] = 'darkorange'
### karyotype colors (http://circos.ca/documentation/tutorials/2d_tracks/connectors/configuration)
map['gpos100'] = rgb(0,0,0, max = 255)
map['gpos'] = rgb(0,0,0, max = 255)
map['gpos75'] = rgb(130,130,130, max = 255)
map['gpos66'] = rgb(160,160,160, max = 255)
map['gpos50'] = rgb(200,200,200, max = 255)
map['gpos33'] = rgb(210,210,210, max = 255)
map['gpos25'] = rgb(200,200,200, max = 255)
map['gvar'] = rgb(220,220,220, max = 255)
map['gneg'] = rgb(255,255,255, max = 255)
map['acen'] = rgb(217,47,39, max = 255)
map['stalk'] = rgb(100,127,164, max = 255)



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

### Structural Variant filter
SV_filter = function(sv_dt, sample_id, min_sv_AD=0, min_sv_span=0) {
    # filter out the SVs w/o any counts in the specified sample
    sv_dt = sv_dt[get(sample_id) > 0]
    # remove 'chr'
    sv_dt[, chr1 := gsub('^chr', '', chr1)]
    sv_dt[, chr2 := gsub('^chr', '', chr2)]
    # generate some new columns
    sv_dt[, `:=`(CONTIG = chr1, POS = pos1, POS2 = pos2), ]
    sv_dt[, SV_color := ifelse((get(sample_id) >= min_sv_AD), 'black', 'gray50')]
    sv_dt[, SV_size := ifelse(chr1 == chr2, pos2 - pos1, -1)]
    sv_dt[, SV_str := str1 * str2]
    sv_dt[, SV_chr2 := ifelse(chr1 != chr2, 0, -1)]
    # create sv1_dt for intrachromosomal SVs
    sv1_dt = sv_dt[chr1 == chr2 & SV_size >= min_sv_span, ]
    sv1_dt = sv1_dt[, .(CONTIG, POS, POS2, SV_color, SV_str, SV_chr2)]
    sv1_dt[, SV_type := 'sv1']
    # create sv2_dt for foldback SVs
    sv2_dt = sv_dt[chr1 == chr2 & SV_size >= min_sv_span & SV_size < 1e6 & SV_str == 1, ]
    sv2_dt = sv2_dt[, .(CONTIG, POS, POS2, SV_color, SV_str, SV_chr2)]
    sv2_dt[, SV_type := 'sv2']
    # create sv3_dt for interchromosomal SVs
    sv31_dt = sv_dt[chr1 != chr2]
    sv32_dt = sv_dt[chr1 != chr2]
    sv31_dt[, `:=`(CONTIG = chr1, POS = pos1, POS2 = pos1, SV_chr2 = chr2), ]
    sv32_dt[, `:=`(CONTIG = chr2, POS = pos2, POS2 = pos2, SV_chr2 = chr1), ]
    sv3_dt = rbind(sv31_dt, sv32_dt)
    sv3_dt = sv3_dt[, .(CONTIG, POS, POS2, SV_color, SV_str, SV_chr2)]
    sv3_dt[, SV_type := 'sv3']
    # combine all the sv_dt
    sv_dt = rbind(sv1_dt, sv2_dt, sv3_dt)
    setkey(sv_dt, CONTIG, POS, POS2)
    # output
    return(sv_dt)
}

### Allele Frequency caller
AF_caller = function(hets_dt, phased_dt) {
    # remove 'chr'
    hets_dt[, CONTIG := gsub('^chr', '', CONTIG)]
    # call AF
    if (nrow(phased_dt) > 0) {
        ## rename colnames for phased VCF
        colnames(phased_dt) = gsub('#', '', colnames(phased_dt))
        colnames(phased_dt)[ncol(phased_dt)] = 'GT'
        phased_dt[, CHROM := gsub('^chr', '', CHROM)]
        phased_dt = phased_dt[, .(CHROM, POS, REF, ALT, GT)]
        ## merge hets TSV and phased VCF
        ephased_dt = merge(hets_dt, phased_dt, all.y = T,
          by.x = c('CONTIG', 'POSITION', 'REF_NUCLEOTIDE', 'ALT_NUCLEOTIDE'), 
          by.y = c('CHROM', 'POS', 'REF', 'ALT'))
        ephased_dt[, c('GT_L', 'GT_R') := tstrsplit(GT, "|", fixed=TRUE)]
        ephased_dt[, XAF := ifelse(GT_L == 1, ALT_COUNT, REF_COUNT)/sum(REF_COUNT, ALT_COUNT), by = .(CONTIG, POSITION)]
        hets_dt = ephased_dt[, .(CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE, XAF)]
    } else {
        hets_dt[, XAF := min(REF_COUNT, ALT_COUNT)/sum(REF_COUNT, ALT_COUNT), by = .(CONTIG, POSITION)]
    }
    # generate a new column
    hets_dt[, POSITION2 := POSITION]
    setkey(hets_dt, CONTIG, POSITION, POSITION2)
    return(hets_dt)
}

### simple Allele Frequency caller
sAF_caller = function(af_dt) {
    saf_dt = copy(af_dt)
    # format
    saf_dt = melt(saf_dt, id=c('CONTIG', 'MIDDLE'), measure=c('XAF1', 'XAF2'), 
                variable.name = 'XAF', value.name = 'MEAN_XAF')
    setkey(saf_dt, CONTIG, MIDDLE)
    return(saf_dt)
}

### Genomic Region extractor
GR_extractor = function(g_region) {
    # remove 'chr'
    g_region = gsub('^chr', '', g_region)
    # split the region
    g_region = unlist(strsplit(g_region, '[:-]'))
    # extract Genomic Region
    if ((length(g_region) == 3) & (as.integer(g_region[2]) <= as.integer(g_region[3]))) {
        g_region_dt = data.table(CONTIG = g_region[1], START_g = as.integer(g_region[2]), END_g = as.integer(g_region[3]))
        setkey(g_region_dt, CONTIG, START_g, END_g)
    } else if (length(g_region) == 1) {
        g_region_dt = data.table(CONTIG = g_region[1])
        setkey(g_region_dt, CONTIG)
    } else {
        cat("ERROR: The value of '-g' must be either 'chr1' or 'chr1:1-100'\n")
        quit(save="no", status=1, runLast=FALSE)
    }
    return(g_region_dt)
}

### Genomic Region clipper
GR_clipper = function(data_dt, ref_dt) {
    if ((ncol(ref_dt) == 3) & (all(ref_dt[['START_g']] <= ref_dt[['END_g']]))) {
        data_dt = foverlaps(data_dt, ref_dt, type = 'within', nomatch=0L)
    } else if (ncol(ref_dt) == 1) {
        data_dt = data_dt[CONTIG %in% ref_dt[['CONTIG']]]
    } else {
        cat("ERROR: The value of '-g' must be either 'chr1' or 'chr1:1-100'\n")
        quit(save="no", status=1, runLast=FALSE)
    }
    return(data_dt)
}

### Genomic Region clipper
GT_extractor = function(data_dt, gene_dt) {
    # process gene table
    colnames(gene_dt) = c('CONTIG', 'START_k', 'ENG_k', 'GENE_NAME', 'GENE_SCORE', 'STRAND_k')
    gene_dt[, CONTIG := gsub('^chr', '', CONTIG)]
    setkey(gene_dt, CONTIG, START_k, ENG_k)
    # overlap genes with copy ratio
    gene_dt = foverlaps(gene_dt, data_dt, type = 'any', nomatch=0L)
    # subset gene table
    gene_dt = gene_dt[, .(CONTIG, START_k, ENG_k, GENE_NAME, GENE_SCORE, STRAND_k)]
    gene_dt = unique(gene_dt)
    # overlap genes with bins
    return(gene_dt)
}

### Copy Number caller
CN_caller = function(cr_dt, hets_dt, gene_dt, purity, ploidy, min_hets_num=0, neu_cn=2) {
    # average ploidy
    D = purity*ploidy + 2*(1-purity)
    # call CN
    if (nrow(hets_dt) > 0) {
        # map the het sites to each bin
        er_dt = foverlaps(hets_dt, cr_dt, type = 'within', nomatch=0L)
        # call mean XAF
        mer_dt = er_dt[, .(MEAN_XAF = mean(XAF, na.rm = TRUE), 
         HETS_N = nrow(.SD[!is.na(XAF)]), 
         ALL_HETS_N = .N), 
          by = .(CONTIG, START, END)]
        cn_dt = merge(cr_dt, mer_dt, by = c('CONTIG', 'START', 'END'))
        cn_dt = cn_dt[MEAN_XAF != 'NaN' & HETS_N >= min_hets_num]
        # call MAJOR_CN and MINOR_CN
        cn_dt[, MAJOR_CN := ((2^LOG2_COPY_RATIO) * D * (1 - MEAN_XAF) - 1) / purity + 1]
        cn_dt[, MINOR_CN := ((2^LOG2_COPY_RATIO) * D * MEAN_XAF - 1) / purity + 1]
        # subset the table
        cn_dt = cn_dt[, .(CONTIG, START, END, MINOR_CN, MAJOR_CN, MEAN_XAF, HETS_N, ALL_HETS_N)]
        #cn_dt = cn_dt[MAJOR_CN > -1 & MINOR_CN > -1]
    } else {
        # call MINOR_CN
        cn_dt = copy(cr_dt)
        cn_dt[, MINOR_CN := (((2^LOG2_COPY_RATIO) * D - 2) / purity + 2) / (2 / 2)]
        cn_dt[CONTIG == 'X', MINOR_CN := (((2^LOG2_COPY_RATIO) * D - 2) / purity + 2) / (2 / neu_cn)]
        # subset the table
        cn_dt = cn_dt[, .(CONTIG, START, END, MINOR_CN)]
        #cn_dt = cn_dt[MINOR_CN > -1]
    }
    # annotate with gene names
    if (nrow(gene_dt) > 0) {
        # call max CN
        mcn_dt = copy(cn_dt)
        mcn_dt[, MAX_CN := apply(.SD, 1, max), .SDcols = grep('_CN$', colnames(mcn_dt), value=TRUE)]
        # overlap genes with bins
        gene_dt = foverlaps(gene_dt, mcn_dt, type = 'any', nomatch=0L)
        # keep the overlaped bins with highest CN
        setorder(gene_dt, -MAX_CN)
        gene_dt = unique(gene_dt, by='GENE_NAME')
        gene_dt = gene_dt[, .(CONTIG, START, END, GENE_NAME, GENE_SCORE)]
        # add gene-related columns
        cn_dt = merge(cn_dt, gene_dt, by = c('CONTIG', 'START', 'END'), all.x = TRUE)
    }
    # determine copy-ratio midpoints
    cn_dt[, MIDDLE := round((START + END) / 2)]
    return(cn_dt)
}

### Allelic Depth caller
AD_caller = function(cr_dt, hets_dt, cn_seg, purity, ploidy, min_hets_num) {
    # average ploidy
    D = purity*ploidy + 2*(1-purity)
    # call TOTAL_CN
    tcn_dt = copy(cr_dt)
    tcn_dt[, TOTAL_CN := ((2^LOG2_COPY_RATIO) * D )]
    tcn_dt = tcn_dt[TOTAL_CN != 'NaN']
    # call TOTAL_COUNT
    tcne_dt = foverlaps(hets_dt, tcn_dt, type = 'within', nomatch=0L)
    tcne_dt[, TOTAL_COUNT := REF_COUNT + ALT_COUNT]
    # call mean values
    tcne_dt = tcne_dt[, .(MEAN_COUNT = mean(TOTAL_COUNT), MEAN_CN = mean(TOTAL_CN), HETS_N = .N),
      by = .(CONTIG, START, END)]
    # filter mean values
    tcne_dt = tcne_dt[MEAN_CN != 'NaN' & HETS_N >= min_hets_num]
    # call AD
    tcne_dt[, AD := MEAN_COUNT / MEAN_CN]
    ad = mean(tcne_dt[['AD']])
    # call threshold of AD
    return(ad)
}

### contig assembler
contig_assembler = function(ref_dict_df, g_region_dt, cn_dt) {
    if (nrow(ref_dict_df) > 0) {
        contig_names = gsub('SN:', '', ref_dict_df[ref_dict_df[,1] == '@SQ', 2])
        contig_lengths = setNames(as.integer(gsub('LN:', '', ref_dict_df[ref_dict_df[,1] == '@SQ', 3])), contig_names)
        if (ncol(g_region_dt) == 3) {
            contig_gl = as.integer(contig_lengths[[g_region_dt[['CONTIG']]]])
            if (g_region_dt[['END_g']] > contig_gl) g_region_dt[['END_g']] = contig_gl
            contig_names = g_region_dt[['CONTIG']]
            contig_ends = g_region_dt[['END_g']] - g_region_dt[['START_g']] + 1
            contig_starts = 0
        } else if (ncol(g_region_dt) == 1) {
            contig_gl = as.integer(contig_lengths[[g_region_dt[['CONTIG']]]])
            contig_names = g_region_dt[['CONTIG']]
            contig_ends = contig_gl
            contig_starts = 0
        } else {
            contig_ends = cumsum(as.numeric(contig_lengths))
            contig_starts = c(0, head(contig_ends, -1))
        }
    } else {
        contig_names = unique(cn_dt[['CONTIG']])
        contig_lengths = cn_dt[, max(END), by = CONTIG][, structure(V1, names = contig_names)]
        if (ncol(g_region_dt) == 3) {
            contig_names = g_region_dt[['CONTIG']]
            contig_ends = g_region_dt[['END_g']] - g_region_dt[['START_g']] + 1
            contig_starts = 0
        } else if (ncol(g_region_dt) == 1) {
            contig_names = g_region_dt[['CONTIG']]
            contig_ends = contig_lengths[g_region_dt[['CONTIG']]]
            contig_starts = 0
        } else {
            contig_ends = cumsum(as.numeric(contig_lengths))
            contig_starts = c(0, head(contig_ends, -1))
        }
    }
    contig_centers = (contig_starts + contig_ends) / 2
    contig_dt = data.table(names = as.factor(contig_names), starts = contig_starts, ends = contig_ends, centers = contig_centers)
    return(contig_dt)
}

### format CN table
CN_format = function(cn_dt, e_dt, p) {
    cn_df = copy(cn_dt)
    if ((nrow(e_dt) > 0) & file.exists(p)) {
        cn_df[, `:=`(H1_CN = MINOR_CN, H2_CN = MAJOR_CN, H1_MEAN_AF = MEAN_XAF), ]
        cn_df = as.data.frame(cn_df[, .(CONTIG, START, END, H1_CN, H2_CN, H1_MEAN_AF, HETS_N, ALL_HETS_N)])
    } else if (nrow(e_dt) > 0) {
        cn_df[, MINOR_MEAN_AF := MEAN_XAF]
        cn_df = as.data.frame(cn_df[, .(CONTIG, START, END, MINOR_CN, MAJOR_CN, MINOR_MEAN_AF)])
    } else {
        cn_df[, TOTAL_CN := MINOR_CN]
        cn_df = as.data.frame(cn_df[, .(CONTIG, START, END, MINOR_CN)])
    }
    return(cn_df)
}


## map
### map cn to contig
cn2contig = function(cn_dt, contig_dt, g_region_dt) {
    # get coordinates
    cn_gdt = copy(cn_dt)
    cn_gdt[, GPOS_X := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + MIDDLE]
    if (ncol(g_region_dt) == 3) cn_gdt[, GPOS_X := GPOS_X - g_region_dt[['START_g']] + 1]
    k_gdt = copy(cn_gdt)
    # transform features to colors
    if ('MAJOR_CN' %in% colnames(cn_gdt)) {
        cn_gdt = melt(cn_gdt, id='GPOS_X', measure=c('MAJOR_CN', 'MINOR_CN'), 
                       variable.name = 'CN_type', value.name = 'CN')
    } else {
        cn_gdt[, CN_type := 'MINOR_CN']
        cn_gdt[, CN := MINOR_CN]
        cn_gdt = cn_gdt[, .(GPOS_X, CN_type, CN)]
    }
    # get coordinates for gene table
    if ('GENE_NAME' %in% colnames(k_gdt)) {
        k_gdt = k_gdt[!is.na(GENE_NAME), .(GPOS_X, GENE_NAME, GENE_SCORE)]
        k_gdt = merge(cn_gdt, k_gdt, by = c('GPOS_X'))
        # keep the overlaped bins with highest CN
        setorder(k_gdt, -CN)
        k_gdt = unique(k_gdt, by='GENE_NAME')
        k_gdt = k_gdt[, .(GPOS_X, CN_type, GENE_NAME, GENE_SCORE)]
        cn_gdt = merge(cn_gdt, k_gdt, by = c('GPOS_X', 'CN_type'), all.x = TRUE)
    }
    cn_gdt[, CN_type := map[unlist(as.character(CN_type))]]
    setorder(cn_gdt, CN_type, GPOS_X)
    return(cn_gdt)
}

### map cn to contig
af2contig = function(cn_dt, contig_dt, g_region_dt) {
    # get coordinates
    af_gdt = copy(cn_dt)
    af_gdt[, GPOS_X := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + MIDDLE]
    if (ncol(g_region_dt) == 3) af_gdt[, GPOS_X := GPOS_X - g_region_dt[['START_g']] + 1]
    # transform features to colors
    af_gdt[, AF_type := 'black']
    af_gdt[, AF := MEAN_XAF]
    af_gdt = af_gdt[, .(GPOS_X, AF_type, AF)]
    return(af_gdt)
}

### map sv to contig
sv2contig = function(sv_dt, contig_dt, g_region_dt) {
    # subset CONTIG
    sv_gdt = sv_dt[CONTIG %in% unique(contig_dt[['names']])]
    if (ncol(g_region_dt) == 3) sv_gdt = foverlaps(sv_gdt, g_region_dt, type = 'within', nomatch=0L)
    if (ncol(g_region_dt) == 1) sv_gdt = sv_gdt[CONTIG == g_region_dt[['CONTIG']]]
    # get coordinates
    sv_gdt[, GPOS_X := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + POS]
    sv_gdt[, GPOS2_X := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + POS2]
    if (ncol(g_region_dt) == 3) {
        sv_gdt[, GPOS_X := GPOS_X - g_region_dt[['START_g']] + 1]
        sv_gdt[, GPOS2_X := GPOS2_X - g_region_dt[['START_g']] + 1]
    }
    sv_gdt = sv_gdt[, .(CONTIG, POS, POS2, SV_type, SV_color, SV_str, SV_chr2, GPOS_X, GPOS2_X)]
    return(sv_gdt)
}

### map sv to contig
cb2contig = function(b_dt, contig_dt, g_region_dt) {
    # subset CONTIG
    b_gdt = b_dt[chr %in% unique(contig_dt[['names']])]
    if ((ncol(g_region_dt) == 1)&(ncol(g_region_dt) == 3)) b_gdt = b_gdt[chr == g_region_dt[['CONTIG']]]
    # get coordinates
    b_gdt[, GPOS_X := contig_dt[['starts']][match(chr, contig_dt[['names']])] + start]
    b_gdt[, GPOS2_X := contig_dt[['starts']][match(chr, contig_dt[['names']])] + end]
    b_gdt[, stain := map[unlist(as.character(stain))]]
    return(b_gdt)
}


## ggplot
### plot SV by geom_bezier (optional)
SV_bezier = function(sv_gdf, contig_df, cn_max) {
    max_xend = max(contig_df[['ends']])
    p = ggplot() + scale_color_identity() + 
      scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01), labels = NULL) + 
      coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
        ylim = c(cn_max, cn_max + 6), expand=FALSE, clip = 'off') + 
      scale_y_continuous(limits = c(cn_max, cn_max + 6), breaks = cn_max + 1.5239, labels = cn_max) + 
      theme(plot.margin = unit(c(2, 2, -3, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = 'bold'), axis.text.y = element_text(size = 10, color = 'white'), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ 
      labs(x = '', y = '')
    sv_gdf[['GPOS1_X']] = (sv_gdf[['GPOS_X']] + sv_gdf[['GPOS2_X']]) / 2
    sv_gdf[['bezier_group']] = paste(sv_gdf[['SV_type']], 1:nrow(sv_gdf), sep = '_')
    sv1_df = sv_gdf[sv_gdf[['SV_type']] == 'sv1', ]
    sv2_df = sv_gdf[sv_gdf[['SV_type']] == 'sv2', ]
    sv3_df = sv_gdf[sv_gdf[['SV_type']] == 'sv3', ]
    if (nrow(sv1_df) > 0) {
        sv1_df = melt(data = sv1_df, 
          id.vars = c('CONTIG', 'SV_color', 'SV_str', 'bezier_group'), 
          measure.vars = c('GPOS_X', 'GPOS1_X', 'GPOS2_X'),
          variable.name = 'GPOS_type', value.name = 'GPOS_X')
        sv1_df = sv1_df[order(sv1_df[['bezier_group']], sv1_df[['GPOS_X']]),]
        sv1_df[['padding_Y']] = 0.1 * (-sv1_df[['SV_str']])
        sv1_df[['GPOS_Y']] = ifelse(sv1_df[['GPOS_type']] == 'GPOS1_X', 2.8, 0) * (-sv1_df[['SV_str']])
        sv1_df[['GPOS_Y']] = sv1_df[['GPOS_Y']] + (cn_max + 3) + sv1_df[['padding_Y']]
        p = p + geom_bezier(aes(x = GPOS_X, y = GPOS_Y, group = bezier_group, color = SV_color), 
          data = sv1_df, size = 0.1, alpha = 1, show.legend = FALSE)
    }
    if (nrow(sv2_df) > 0) {
        sv2_df = melt(data = sv2_df, 
          id.vars = c('CONTIG', 'SV_color', 'SV_str', 'bezier_group'), 
          measure.vars = c('GPOS_X', 'GPOS1_X', 'GPOS2_X'),
          variable.name = 'GPOS_type', value.name = 'GPOS_X')
        sv2_df = sv2_df[order(sv2_df[['bezier_group']], sv2_df[['GPOS_X']]),]
        sv2_df[['padding_Y']] = 0.1 * (-sv2_df[['SV_str']])
        sv2_df[['GPOS_Y']] = ifelse(sv2_df[['GPOS_type']] == 'GPOS1_X', 2.8, 0) * (-sv2_df[['SV_str']])
        sv2_df[['GPOS_Y']] = sv2_df[['GPOS_Y']] + (cn_max + 3) + sv2_df[['padding_Y']]
        sv2p_df = sv2_df[(sv2_df[['SV_str']] == 1)&(sv2_df[['GPOS_type']] != 'GPOS1_X'),]
        if (nrow(sv2p_df) > 0) {
            sv2p_df[['SV_color']] = ifelse(sv2p_df[['SV_color']] == 'black', 'violet', 'pink')
            p = p + geom_segment(aes(x = GPOS_X, y = cn_max + 2.9, xend = GPOS_X, yend = cn_max + 1.5, color = SV_color), 
                  data = sv2p_df, size = 0.1, alpha = 1, show.legend = FALSE)
        }
    }
    if (nrow(sv3_df) > 0) {
        sv3_df[['SV_color']] = ifelse(sv3_df[['SV_color']] == 'black', 'violet', 'pink')
        p = p + geom_segment(aes(x = GPOS_X, y = cn_max + 3.1, xend = GPOS2_X, yend = cn_max + 4.5, color = SV_color), 
          data = sv3_df, size = 0.1, alpha = 1, show.legend = FALSE)
        p = p + geom_text_repel(aes(x = GPOS_X, y = cn_max + 4.5, label = SV_chr2, color = SV_color), 
          data = sv3_df, segment.size = 0.1, force = 0.1, 
          direction = 'x', nudge_x = 1, nudge_y = 0.5, box.padding = 0, show.legend = FALSE)
    }
    return(p)
}

### plot CN by geom_point
CN_point = function(cn_gdf, contig_df, cn_max) {
    max_xend = contig_df[nrow(contig_df), 'ends']
    contig_df[nrow(contig_df), 'ends'] = max_xend + max_xend * .01
    cn_gdf[cn_gdf[['CN']] > 2, 'CN'] = log2(cn_gdf[cn_gdf[['CN']] > 2, 'CN']) + 1
    p = ggplot() + scale_color_identity() + 
      geom_point(aes(x = GPOS_X, y = CN, color = CN_type), 
        data = cn_gdf, size = 0.5, show.legend = FALSE) + 
      geom_hline(yintercept=c(0, 0.5, 1, 1.5, 2, log2(c(3, 4)) + 1), linetype='dashed', color = 'gray') + 
      geom_segment(aes(x = ends, y = -1, xend = ends, yend = log2(cn_max) + 1), 
        data = contig_df, color = 'gray80', show.legend = FALSE) + 
      scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01)) + 
      coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
        ylim = c(-1, log2(cn_max) + 1), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(-1, log2(cn_max) + 1), 
        breaks = seq(-1, log2(cn_max) + 1), labels = c(-1, 0, 2 ^ seq(0, log2(cn_max)))) + 
      theme(plot.margin = unit(c(0, 2, 0, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), axis.line.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = 'bold'), axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    if (length(unique(cn_gdf[['CN_type']])) == 2) p = p + labs(x = '', y = 'Allelic CN')
    if (length(unique(cn_gdf[['CN_type']])) == 1) p = p + labs(x = '', y = 'Total CN')
    if ('GENE_NAME' %in% colnames(cn_gdf)) {
        p = p + 
          geom_text_repel(aes(x = GPOS_X, y = CN, label = GENE_NAME, color = GENE_SCORE), 
            data = cn_gdf[!is.na(cn_gdf[['GENE_NAME']]), ], segment.color = 'black', segment.alpha = 0.5,
            nudge_x = 1, nudge_y = 1, show.legend = FALSE)
    }
    return(p)
}

### plot AF by geom_point
AF_point = function(af_gdf, contig_df, af_max) {
    max_xend = contig_df[nrow(contig_df), 'ends']
    contig_df[nrow(contig_df), 'ends'] = max_xend + max_xend * .01
    p = ggplot() + scale_color_identity() + 
      geom_point(aes(x = GPOS_X, y = AF, color = AF_type), 
        data = af_gdf, size = 0.5, show.legend = FALSE) + 
      geom_hline(yintercept=c(0, 0.5, 1), linetype='dashed', color = 'gray') + 
      geom_segment(aes(x = ends, y = 0, xend = ends, yend = af_max), 
        data = contig_df, color = 'gray80', show.legend = FALSE) + 
      scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01)) + 
      coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
        ylim = c(0, af_max), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(0, af_max), 
        breaks = seq(0, af_max), labels = c(0, af_max)) + 
      theme(plot.margin = unit(c(0, 2, 0, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), axis.line.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = 'bold'), axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
        labs(x = '', y = 'AF')
    if ('GENE_NAME' %in% colnames(af_gdf)) {
        p = p + 
          geom_text_repel(aes(x = GPOS_X, y = AF, label = GENE_NAME, color = GENE_SCORE), 
            data = af_gdf[!is.na(af_gdf[['GENE_NAME']]), ], segment.color = 'black', segment.alpha = 0.5,
            nudge_x = 1, nudge_y = 1, show.legend = FALSE)
    }
    return(p)
}

### plot cytoband by geom_rect
cytoband_rect = function(b_gdf, contig_df, g_region_dt, sn) {
    p = ggplot() + scale_fill_identity() + 
      scale_y_continuous(limits = c(-2, -1)) + 
      theme(plot.margin = unit(c(0, 2, 2, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), axis.line.y = element_blank(), 
        axis.title.x = element_text(size = 12, face = 'bold'), axis.text.x = element_text(size = 10), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    if ((nrow(contig_df) == 1)&(ncol(g_region_dt) == 1)) {
        max_xend = contig_df[nrow(contig_df), 'ends']
        init_x = as.integer(substr(max_xend, 1, 1))
        scale_x = ifelse(init_x > 2, 10^(floor(log10(max_xend))), 10^(floor(log10(max_xend)) - 1))
        unit_x = unit_converter(scale_x)
        unit_x = substring(unit_x, nchar(unit_x), nchar(unit_x))
        breaks_x = seq(0, max_xend, scale_x)
        p = p + 
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -1.67, ymax = -1.25, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.5, show.legend = FALSE) + 
          coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            ylim = c(-2, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            breaks = breaks_x, 
            labels = gsub(paste0(' ', unit_x), '', unit_converter(breaks_x))) + 
          labs(x = paste('Chromosome', g_region_dt[['CONTIG']][1], paste0('(', unit_x, 'b)'), 'in', sn), y = '')
    } else if ((nrow(contig_df) == 1)&(ncol(g_region_dt) == 3)) {
        max_xend = b_gdf[nrow(b_gdf), 'end']
        init_x = as.integer(substr(max_xend, 1, 1))
        scale_x = ifelse(init_x > 2, 10^(floor(log10(max_xend))), 10^(floor(log10(max_xend)) - 1))
        unit_x = unit_converter(scale_x)
        unit_x = substring(unit_x, nchar(unit_x), nchar(unit_x))
        breaks_x = seq(0, max_xend, scale_x)
        p = p + 
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -1.67, ymax = -1.25, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.5, show.legend = FALSE) + 
          coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            ylim = c(-2, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            breaks = breaks_x, 
            labels = gsub(paste0(' ', unit_x), '', unit_converter(breaks_x))) + 
          geom_rect(aes(xmin = START_g, xmax = END_g, ymin = -1.65, ymax = -1.25), 
            data = g_region_dt, fill = 'red', colour = 'red', alpha = 0.5, size = 0.5, show.legend = FALSE) + 
          labs(x = paste('Chromosome', paste0(g_region_dt[['CONTIG']][1], ':', g_region_dt[['START_g']][1], '-', g_region_dt[['END_g']][1]), paste0('(', unit_x, 'b)'), 'in', sn), y = '')
    } else {
        max_xend = contig_df[nrow(contig_df), 'ends']
        p = p + 
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -1.67, ymax = -1.25, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.1, show.legend = FALSE) + 
          geom_segment(aes(x = ends, y = -2, xend = ends, yend = -1), 
            data = head(contig_df, -1), color = 'gray80', show.legend = FALSE) + 
          coord_cartesian(xlim = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            ylim = c(-2, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0 - max_xend * .01, max_xend + max_xend * .01), 
            breaks = contig_df[['centers']], labels = contig_df[['names']]) + 
          labs(x = paste('Chromosomes', 'in', sn), y = '')
    }
    return(p)
}

### generate bottom panel of the ouput plot
CN_contour = function(cn_df, sv_df = NULL) {
    cn4 = colnames(cn_df)[4]
    cn5 = colnames(cn_df)[5]
    cn4_max = 2^ceiling(log2(max(cn_df[, cn4])))
    cn5_max = 2^ceiling(log2(max(cn_df[, cn5])))
    cn_df[cn_df[[cn4]] > 2, cn4] = log2(cn_df[cn_df[[cn4]] > 2, cn4]) + 1
    cn_df[cn_df[[cn5]] > 2, cn5] = log2(cn_df[cn_df[[cn5]] > 2, cn5]) + 1
    cn2d = MASS::kde2d(cn_df[[cn4]], cn_df[[cn5]], n = 100)
    cne_df = expand.grid(cn4e = cn2d$x, cn5e = cn2d$y)
    cne_df[['dens']] = as.vector(cn2d$z)
    cne_bin = max(as.vector(cn2d$z)) * c(.99, .95, .9, .75, .5, .25, .1)
    p = ggplot() + scale_color_continuous(low = 'blue', high = 'cyan') + 
      geom_point(aes(x = get(cn4), y = get(cn5)), 
        data = cn_df, color = 'black', alpha = 0.5, size = 0.5, show.legend = FALSE) + 
      geom_contour(aes(x = cn4e, y = cn5e, z = dens, colour = ..level..), 
        data = cne_df, breaks = cne_bin, size= 0.5) + 
      geom_vline(xintercept=c(0, 1, 2, log2(c(3, 4)) + 1), linetype='dashed', color = 'gray') + 
      geom_hline(yintercept=c(0, 1, 2, log2(c(3, 4)) + 1), linetype='dashed', color = 'gray') + 
      coord_cartesian(xlim = c(-1, log2(cn4_max) + 1), 
        ylim = c(-1, log2(cn5_max) + 1), expand = FALSE, clip = 'off') +
      scale_x_continuous(limits = c(-1, log2(cn4_max) + 1), 
        breaks = seq(-1, log2(cn4_max) + 1), labels = c(-1, 0, 2 ^ seq(0, log2(cn4_max)))) + 
      scale_y_continuous(limits = c(-1, log2(cn5_max) + 1), 
        breaks = seq(-1, log2(cn5_max) + 1), labels = c(-1, 0, 2 ^ seq(0, log2(cn5_max)))) + 
      theme(plot.margin = unit(c(2, 2, 2, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), legend.position = 'none', 
        axis.title = element_text(size = 12, face = 'bold'), axis.text = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
      labs(x = gsub('_', ' ', cn4), y = gsub('_', ' ', cn5))
    return(p)
}

### generate layout of CN plots
CN_scatter = function(cn_gdf, sv_gdf, contig_df) {
    cn_max = 2^ceiling(log2(max(cn_gdf[['CN']])))
    if (nrow(sv_gdf) > 0) {
        sva_p = SV_bezier(sv_gdf=sv_gdf, contig_df=contig_df, cn_max=cn_max)
        cnl_p = CN_point(cn_gdf=cn_gdf, contig_df=contig_df, cn_max=cn_max)
        g0 <- ggplotGrob(sva_p); g2 <- ggplotGrob(cnl_p)
        g0$heights[unique(g0$layout[g0$layout$name == 'panel', 't'])] = unit(6, 'cm')
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(3 + (ceiling(log2(cn_max))-1)/2*2, 'cm')
        gb <- rbind(g0, g2, size='first')
    } else {
        cnl_p = CN_point(cn_gdf=cn_gdf, contig_df=contig_df, cn_max=cn_max)
        g2 <- ggplotGrob(cnl_p)
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(3 + (ceiling(log2(cn_max))-1)/2*2, 'cm')
        gb <- rbind(g2, size='first')
    }
    if (nrow(b_gdf) > 0) {
        ground_p = cytoband_rect(b_gdf=b_gdf, contig_df=contig_df, g_region_dt=g_df, sn=s)
        g3 <- ggplotGrob(ground_p)
        g3$heights[unique(g3$layout[g3$layout$name == 'panel', 't'])] = unit(1, 'cm')
        gb <- rbind(gb, g3, size='first')
    }
    return(gb)
}

### generate layout of AF plots
AF_scatter = function(af_gdt, sv_gdf, contig_df) {
    af_max = 1
    if (nrow(sv_gdf) > 0) {
        sva_p = SV_bezier(sv_gdf=sv_gdf, contig_df=contig_df, cn_max=af_max)
        afl_p = AF_point(af_gdf=af_gdt, contig_df=contig_df, af_max=af_max)
        g0 <- ggplotGrob(sva_p); g2 <- ggplotGrob(afl_p)
        g0$heights[unique(g0$layout[g0$layout$name == 'panel', 't'])] = unit(6, 'cm')
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(6, 'cm')
        gb <- rbind(g0, g2, size = 'first')
    } else {
        afl_p = AF_point(af_gdf=af_gdt, contig_df=contig_df, af_max=af_max)
        g2 <- ggplotGrob(afl_p)
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(6, 'cm')
        gb <- rbind(g2, size = 'first')
    }
    if (nrow(b_gdf) > 0) {
        ground_p = cytoband_rect(b_gdf = b_gdf, contig_df = contig_df, g_region_dt = g_df, sn = s)
        g3 <- ggplotGrob(ground_p)
        g3$heights[unique(g3$layout[g3$layout$name == 'panel', 't'])] = unit(1, 'cm')
        gb <- rbind(gb, g3, size = 'first')
    }
    return(gb)
}

### generate layout of AF plots
CN_scatter2 = function(cn_df, xmax, ymax) {
    cnh_p = CN_contour(cn_df = cn_df)
    gb <- ggplotGrob(cnh_p)
    gb$widths[unique(gb$layout[gb$layout$name == 'panel', 'r'])] = unit(3 + (ceiling(log2(xmax))-1)/2*2, 'cm')
    gb$heigbts[unique(gb$layout[gb$layout$name == 'panel', 't'])] = unit(3 + (ceiling(log2(ymax))-1)/2*2, 'cm')
    return(gb)
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
    make_option(c("-s", "--sample_name"), type = "character", default = "tumor",
     help = "Sample name for output [%default]", metavar = "character"),
    make_option(c("-d", "--ref_dict"), type = "character", default = "",
     help = "Reference file in DICT format [%default]", metavar = "character"),
    make_option(c("-k", "--gene_bed"), type = "character", default = "",
     help = "Gene file in BED format [%default]", metavar = "character"),
    make_option(c("-n", "--sv_tsv"), type = "character", default = "",
     help = "Structural variation data in bkps-TSV format. (optinal) [%default]", metavar = "character"),
    make_option(c("-e", "--hets_tsv"), type = "character", default = "",
     help = "Het sites data file in TSV format. The first four columns must be \n\t\tCONTIG, POSITION, REF_COUNT and ALT_COUNT. (optinal) [%default]", metavar = "character"),
    make_option(c("-p", "--phased_vcf"), type = "character", default = "",
     help = "Phased genotypes in VCF format (single sample VCF). Please refer to 'https://imputation.sanger.ac.uk/' (optinal) [%default]", metavar = "character"),
    make_option(c("-b", "--cytoband_tsv"), type = "character", default = "",
     help = "Cytoband data in TSV format. Please refer to 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz' (optinal) [%default]", metavar = "character"),
    make_option(c("-u", "--purity"), type = "double", default = "1",
     help = "Purity [%default]", metavar = "double"),
    make_option(c("-l", "--ploidy"), type = "double", default = "2",
     help = "Ploidy [%default]", metavar = "double"),
    make_option(c("-g", "--region"), type = "character", default = "",
     help = "The output region [%default] \n\t\te.g. g='chr17:36000000-42000000'", metavar = "character"),
    make_option(c("-x", "--neu_cn"), type = "integer", default = "2",
     help = "Copy number of X chromosome [%default] \n\t\te.g. 2 for female's X, 1 for male's X", metavar = "integer"),
    make_option(c("--th_ad"), type = "double", default = "-1",
     help = "Threshold for allelic depth in the specified sample [%default] \n\t\tIf 'hets_tsv' exists, '-1' means this threshold will be infered from the data. \n\t\tIf 'hets_tsv' does not exist, '-1' will be converted to 3.", metavar = "double"),
    make_option(c("--th_span"), type = "double", default = "1e5",
     help = "Threshold for SV span [%default]", metavar = "double"),
    make_option(c("--th_hetn"), type = "double", default = "-1",
     help = "Threshold for hets number (roughly, 0 for 5kb, 2 for 10kb, 5 for 25kb) [%default]", metavar = "double"),
    make_option(c("--plot_width"), type = "double", default = "12",
     help = "Width of the ouput plot [%default]", metavar = "double"),
    make_option(c("--plot_height"), type = "double", default = "7",
     help = "Height of the ouput plot [%default]", metavar = "double"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "plotACN",
     help = "Output prefix [%default]", metavar = "character"))

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_ratio.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The first four columns of <copy_ratio.tsv> must be CONTIG, START, END and LOG2_COPY_RATIO.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

s = opt$options$sample_name
d = opt$options$ref_dict
k = opt$options$gene_bed
n = opt$options$sv_tsv
e = opt$options$hets_tsv
p = opt$options$phased_vcf
b = opt$options$cytoband_tsv
u = opt$options$purity
l = opt$options$ploidy
g = opt$options$region
x = opt$options$neu_cn
pw = opt$options$plot_width
ph = opt$options$plot_height
o = opt$options$output_prefix
ta = opt$options$th_ad
tp = opt$options$th_span
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
    n = ifelse(n == '-', 'file:///dev/stdin', n)
    p = ifelse(p == '-', 'file:///dev/stdin', p)
    if (file.exists(p)) {if(grepl('\\.gz$', p)) p_dt = fread(cmd = paste0('zcat ', p), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
    if (file.exists(b)) {if(grepl('\\.gz$', b)) b_dt = fread(cmd = paste0('zcat ', b), skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
} else {
    r = ifelse(r == '-', 'cat /dev/stdin', r)
    e = ifelse(e == '-', 'cat /dev/stdin', e)
    n = ifelse(n == '-', 'cat /dev/stdin', n)
    p = ifelse(p == '-', 'cat /dev/stdin', p)
    if (file.exists(p)) {if(grepl('\\.gz$', p)) p_dt = fread(paste0('zcat ', p), skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
    if (file.exists(b)) {if(grepl('\\.gz$', b)) b_dt = fread(paste0('zcat ', b), skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
}
### fread
r_dt = fread(r, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(e)) e_dt = fread(e, skip = 'CONTIG', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(n)) sv_dt = fread(n, skip = 'SVidx', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(p)) {if(!grepl('\\.gz$', p)) p_dt = fread(p, skip = '#CHROM', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
if (file.exists(b)) {if(!grepl('\\.gz$', b)) b_dt = fread(b, skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
if (file.exists(d)) contig_dict_df = read.delim(d, header=FALSE)
if (file_non0(k)) gene_dt = fread(k, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = FALSE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)


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
