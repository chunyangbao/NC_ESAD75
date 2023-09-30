#!/usr/bin/env Rscript

# functions
## color
### color vector
map = NULL
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

### Copy Number extractor
CN_extractor = function(cn_dt) {
    if ('PH2_CN' %in% colnames(cn_dt)) {
        if(nrow(cn_dt[!is.numeric(PH2_CN)]) > nrow(cn_dt)/2) {
            cn_dt = cn_dt[, .(CONTIG, START, END, PH1_CN)]
        }else{
        cn_dt = cn_dt[, .(CONTIG, START, END, PH1_CN, PH2_CN)]
        }
    } else {
        colnames(cn_dt)[grep('CN', colnames(cn_dt))[1]] = 'PH1_CN'
        cn_dt = cn_dt[, .(CONTIG, START, END, PH1_CN)]
    }
    setkey(cn_dt, CONTIG, START, END)
    return(cn_dt)
}

### Structural Variation extractor
SV_extractor = function(sv_dt) {
    sv_dt = sv_dt[, .(CONTIG, POS, POS2, SV_type, SV_color, SV_str, SV_chr2)]
    setkey(sv_dt, CONTIG, POS, POS2)
    return(sv_dt)
}

### Cytoband extractor
CB_extractor = function(b_dt) {
    b_dt[, chr := gsub('^chr', '', chr)]
    setkey(b_dt, chr, start, end)
    return(b_dt)
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

### segments extractor
SG_extractor = function(seg_dt) {
    seg_dt[, SAMPLE := gsub('^.*\\.PH', 'PH', SAMPLE)]
    seg_dt[, SAMPLE := gsub('\\..*$', '', SAMPLE)]
    seg_dt[, CONTIG := as.character(CONTIG)]
    setkey(seg_dt, CONTIG, START, END)
    return(seg_dt)
}

### Genomic Region clipper
SG_clipper = function(data_dt, ref_dt) {
    if ((ncol(ref_dt) == 3) & (all(ref_dt[['START_g']] <= ref_dt[['END_g']]))) {
        ref_start = as.integer(ref_dt[1, 'START_g'])
        ref_end = as.integer(ref_dt[1, 'END_g'])
        data_dt = foverlaps(data_dt, ref_dt, type = 'any', nomatch=0L)
        data_dt[START <= ref_start, START := ref_start + 1]
        data_dt[END >= ref_end, END := ref_end]
    } else if (ncol(ref_dt) == 1) {
        data_dt = data_dt[CONTIG %in% ref_dt[['CONTIG']]]
    } else {
        cat("ERROR: The value of '-g' must be either 'chr1' or 'chr1:1-100'\n")
        quit(save="no", status=1, runLast=FALSE)
    }
    data_dt = data_dt[, .(SAMPLE, CONTIG, START, END, BINS, CN)]
    setkey(data_dt, SAMPLE, CONTIG, START, END)
    return(data_dt)
}

CN_caller = function(cn_dt, gene_dt) {
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


## map
### map cn to contig
cn2contig = function(cn_dt, contig_dt, g_region_dt) {
    # get coordinates
    cn_gdt = copy(cn_dt)
    cn_gdt[, GPOS_X := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + MIDDLE]
    if (ncol(g_region_dt) == 3) cn_gdt[, GPOS_X := GPOS_X - g_region_dt[['START_g']] + 1]
    k_gdt = copy(cn_gdt)
    # transform features to colors
    if ('PH2_CN' %in% colnames(cn_gdt)) {
        cn_gdt = melt(cn_gdt, id=c('CONTIG', 'START', 'END', 'GPOS_X'), measure=c('PH1_CN', 'PH2_CN'), 
                       variable.name = 'CN_type', value.name = 'CN')
    } else {
        cn_gdt[, CN_type := 'PH1_CN']
        cn_gdt[, CN := PH1_CN]
        cn_gdt = cn_gdt[, .(CONTIG, START, END, GPOS_X, CN_type, CN)]
    }
    # get coordinates for gene table
    if ('GENE_NAME' %in% colnames(k_gdt)) {
        k_gdt = k_gdt[!is.na(GENE_NAME), .(GPOS_X, GENE_NAME, GENE_SCORE)]
        k_gdt = merge(cn_gdt, k_gdt, by = c('GPOS_X'))
        # keep the overlaped bins with highest CN
        setorder(k_gdt, -CN)
        k_gdt = unique(k_gdt, by='GENE_NAME')
        k_gdt = k_gdt[, .(GPOS_X, CN_type, GENE_NAME, GENE_SCORE)]
        cn_gdt = merge(cn_gdt, k_gdt, by = c('GPOS_X', 'CN_type'), all.x = TRUE, allow.cartesian = TRUE)
        cn_gdt = unique(cn_gdt, by=c('GENE_NAME', 'CN_type', 'GPOS_X'))
    }
    cn_gdt[, CN_color := gsub('CN', 'color', CN_type)]
    cn_gdt[, CN_shape := gsub('CN', 'shape', CN_type)]
    cn_gdt[, CN_color := map[unlist(as.character(CN_color))]]
    cn_gdt[, CN_shape := map[unlist(as.character(CN_shape))]]
    setorder(cn_gdt, -CN_type, CONTIG, GPOS_X)
    return(cn_gdt)
}

### map cn_seg to contig
cns2contig = function(cn_seg, contig_dt, g_region_dt) {
    # get coordinates
    cn_gseg = copy(cn_seg)
    cn_gseg[, GPOS_X0 := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + START]
    cn_gseg[, GPOS_X1 := contig_dt[['starts']][match(CONTIG, contig_dt[['names']])] + END]
    if (ncol(g_region_dt) == 3) cn_gseg[, GPOS_X0 := GPOS_X0 - g_region_dt[['START_g']] + 1]
    if (ncol(g_region_dt) == 3) cn_gseg[, GPOS_X1 := GPOS_X1 - g_region_dt[['START_g']] + 1]
    cn_gseg[, SAMPLE := map[unlist(as.character(SAMPLE))]]
    return(cn_gseg)
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
      scale_x_continuous(limits = c(0, max_xend + max_xend * .01), labels = NULL) + 
      coord_cartesian(xlim = c(0, max_xend), 
        ylim = c(cn_max, cn_max + 6), expand=FALSE, clip = 'off') + 
      scale_y_continuous(limits = c(cn_max, cn_max + 6)) + 
      theme(plot.margin = unit(c(0, 0, -1.5, 0), 'cm'), plot.background = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
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

### call change points for CN
CN_change = function(cn_seg_gdf, cytoband_gdf) {
    cn_seg_gdft = data.table(cn_seg_gdf)
    # annotate with SAMPLE
    s1 = cn_seg_gdft[1, SAMPLE]
    cn_seg_gdft[, is_s1 := ifelse(SAMPLE == s1, 1, 0)]
    # call changing points
    cn_seg_gdft[, s1_cp := c(diff(is_s1), 1)]
    cn_seg_gdft[, diff_CN := c(diff(CN), 0)]
    cn_seg_gdft[, CN2 := CN + diff_CN]
    setkey(cn_seg_gdft, CONTIG, START, END)
    # remove the one closed to acen
    if (nrow(cytoband_gdf) > 0) {
        cytoband_gdt = data.table(cytoband_gdf)
        acen_gdt = cytoband_gdt[stain == '#D92F27']
        # find targets (assumption: each chr has one acen)
        acen_gdt[, `:=`(CONTIG = chr, START_a = 0, END_a = start)]
        acen_gdt = acen_gdt[, .(CONTIG, START_a, END_a)]
        setkey(acen_gdt, CONTIG, START_a, END_a) 
        acen_gdt = unique(acen_gdt, by='CONTIG') # determin the region of p-arm
        acen_gdt = foverlaps(cn_seg_gdft, acen_gdt, nomatch=NULL)
        if(nrow(acen_gdt) > 0) {
            acen_gdt = acen_gdt[, max(END), by=.(SAMPLE, CONTIG)]
        } else {
            acen_gdt = acen_gdt[, sum(END), by=.(SAMPLE, CONTIG)]
        }
        acen_gdt[, `:=`(END = V1, is_acen = 1)]
        acen_gdt = acen_gdt[, .(SAMPLE, CONTIG, END, is_acen)]
        # remove the targets
        cn_seg_gdft = merge(cn_seg_gdft, acen_gdt, by = c('SAMPLE', 'CONTIG', 'END'), all.x = TRUE)
        cn_seg_gdft = cn_seg_gdft[is.na(is_acen)]
        cn_seg_gdft[, is_acen := NULL]
    }
    # output
    cn_seg_gdft = cn_seg_gdft[(s1_cp == 0), .(SAMPLE, CONTIG, START, END, GPOS_X1, CN, CN2)]
    setDF(cn_seg_gdft)
    return(cn_seg_gdft)
}

### plot CN by geom_point
#https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
#https://github.com/r-lib/gtable
CN_point = function(cn_gdf, cn_seg_gdf, contig_df, cytoband_gdf, cn_max, dots_size) {
    max_xend = contig_df[nrow(contig_df), 'ends']
    cn_gdf[cn_gdf[['CN']] > 2, 'CN'] = log2(cn_gdf[cn_gdf[['CN']] > 2, 'CN']) + 1
    p = ggplot() + scale_color_identity() + scale_fill_identity() + scale_shape_identity() + 
      geom_point(aes(x = GPOS_X, y = CN, color = CN_color, fill = CN_color, shape = CN_shape), 
        data = cn_gdf, size = dots_size, alpha = 1, show.legend = FALSE) + 
      geom_hline(yintercept=c(0, 0.5, 1, 1.5, 2, log2(c(3, 4)) + 1), linetype='dashed', color = 'black', size = 0.4) + 
      scale_x_continuous(limits = c(0, max_xend + max_xend * .01)) + 
      coord_cartesian(xlim = c(0, max_xend), 
        ylim = c(0, log2(cn_max) + 1), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(-1, log2(cn_max) + 1), 
        breaks = seq(0, log2(cn_max) + 1), labels = c(0, 2 ^ seq(0, log2(cn_max)))) + 
      theme(plot.margin = unit(c(0, 1, 1, 0), 'cm'), plot.background = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    # add segments
    if (nrow(cn_seg_gdf) > 0) {
        cn_seg_gdf[cn_seg_gdf[['CN']] > 2, 'CN'] = log2(cn_seg_gdf[cn_seg_gdf[['CN']] > 2, 'CN']) + 1
        cnc_gdf = CN_change(cn_seg_gdf, cytoband_gdf)
        p = p + 
          geom_segment(aes(x = GPOS_X0, y = CN, xend = GPOS_X1, yend = CN, color = SAMPLE), 
            data = cn_seg_gdf, alpha = 0.5, size = 1, show.legend = FALSE) + 
          geom_segment(aes(x = GPOS_X1 + 0.5, y = CN, xend = GPOS_X1 + 0.5, yend = CN2, color = SAMPLE), 
            data = cnc_gdf, alpha = 0.5, size = 0.4, show.legend = FALSE)
    }
    # add gene names
    if ('GENE_NAME' %in% colnames(cn_gdf)) {
        p = p + 
          geom_text_repel(aes(x = GPOS_X, y = CN, label = GENE_NAME, color = GENE_SCORE), 
            data = cn_gdf[!is.na(cn_gdf[['GENE_NAME']]), ], segment.color = 'black', segment.alpha = 0.5, size = 3,
            nudge_x = 1, nudge_y = 1, show.legend = FALSE)
    }
    # add a line by the end of each chr for multi-chr plot
    if (nrow(contig_df) > 1) {
        contig_df[nrow(contig_df), 'ends'] = max_xend + max_xend * .01
        p = p +
          geom_segment(aes(x = ends, y = -1, xend = ends, yend = log2(cn_max) + 1), 
            data = contig_df, color = 'gray80', show.legend = FALSE)
    }

    return(p)
}

### plot x-axis
xaxis_bottom = function(contig_df, region_df) {
    # get the length
    xstart = region_df[1, 2]
    max_xend = contig_df[nrow(contig_df), 'ends']
    # get the scale
    scale_x = 10^(floor(log10(region_df[1, 3])))
    unit_x = unit_converter(scale_x)
    unit_x = substring(unit_x, nchar(unit_x), nchar(unit_x))
    breaks_x = seq(region_df[1, 2], region_df[1, 3], scale_x)
    # plot
    p = ggplot() + 
      geom_point(aes(x = 0, y = 0), alpha = 0, show.legend = FALSE) + 
      scale_x_continuous(limits = c(0, max_xend + max_xend * .01), 
        breaks = breaks_x - xstart, labels = gsub(paste0(' ', unit_x), '', unit_converter(breaks_x))) + 
      coord_cartesian(xlim = c(0, max_xend), 
        ylim = c(0, 0), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(0, 0)) + 
      theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'), plot.background = element_blank(), 
        axis.ticks.length = unit(0.3, 'cm'),
        axis.line.x = element_line(colour = 'black', size = 0.5), axis.line.y = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_text(size = 18, colour = 'black'), 
        axis.ticks.x = element_line(colour = 'black', size = 0.5), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    return(p)
}

### plot y-axis for CN_point
yaxis_left = function(cn_gdf, cn_max, ylab) {
    p = ggplot() + 
      geom_point(aes(x = 0, y = 0), alpha = 0, show.legend = FALSE) + 
      scale_x_continuous(limits = c(0, 0), labels = NULL) + 
      coord_cartesian(xlim = c(0, 0), 
        ylim = c(0, log2(cn_max) + 1), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(-1, log2(cn_max) + 1), 
        breaks = seq(0, log2(cn_max) + 1), labels = c(0, 2 ^ seq(0, log2(cn_max)))) + 
      theme(plot.margin = unit(c(0, 0, 1, 1), 'cm'), plot.background = element_blank(), 
        axis.ticks.length = unit(0.3, 'cm'),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.title.y = element_text(size = 18, face = 'bold', colour = 'black'), 
        axis.text.y = element_text(size = 18, colour = 'black'), 
        axis.ticks.y = element_line(colour = 'black', size = 0.5), 
        axis.line = element_line(colour = 'black', size = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    p = p + labs(x = '', y = ylab)
    return(p)
}

### plot y-axis for CN_point
yaxis0_left = function(ymin, ymax) {
    p = ggplot() + 
      scale_x_continuous(limits = c(0, 0), labels = NULL) + 
      coord_cartesian(xlim = c(0, 0), 
        ylim = c(ymin, ymax), expand = FALSE, clip = 'off') +
      scale_y_continuous(limits = c(ymin, ymax)) + 
      theme(plot.margin = unit(c(0, 0, 0, 1), 'cm'), plot.background = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    if (length(unique(cn_gdf[['CN_type']])) == 2) p = p + labs(x = '', y = 'X')
    return(p)
}

### plot cytoband by geom_rect
cytoband_rect = function(cytoband_gdf, contig_df, g_region_dt) {
    b_gdf = cytoband_gdf
    p = ggplot() + scale_fill_identity() + 
      scale_y_continuous(limits = c(-2, -1)) + 
      theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'), plot.background = element_blank(), 
        axis.ticks.length = unit(0.3, 'cm'),
        axis.line = element_line(colour = 'black', size = 0.5), axis.line.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 18, colour = 'black'), 
        axis.ticks.x = element_line(colour = 'black', size = 0.5), 
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
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -2, ymax = -1.4, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.5, show.legend = FALSE) + 
          coord_cartesian(xlim = c(0, max_xend), 
            ylim = c(-2.4, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0, max_xend + max_xend * .01), 
            breaks = breaks_x, 
            labels = gsub(paste0(' ', unit_x), '', unit_converter(breaks_x)))
    } else if ((nrow(contig_df) == 1)&(ncol(g_region_dt) == 3)) {
        max_xend = b_gdf[nrow(b_gdf), 'end']
        init_x = as.integer(substr(max_xend, 1, 1))
        scale_x = ifelse(init_x > 2, 10^(floor(log10(max_xend))), 10^(floor(log10(max_xend)) - 1))
        unit_x = unit_converter(scale_x)
        unit_x = substring(unit_x, nchar(unit_x), nchar(unit_x))
        breaks_x = seq(0, max_xend, scale_x)
        p = p + 
         # geom_point(aes(x = c(0 - max_xend * .01, 0 + max_xend * .01), y = 0), alpha = 1, show.legend = FALSE) + 
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -2, ymax = -1.4, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.5, show.legend = FALSE) + 
          coord_cartesian(xlim = c(0, max_xend), 
            ylim = c(-2.4, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0, max_xend + max_xend * .01), 
            breaks = breaks_x, 
            labels = gsub(paste0(' ', unit_x), '', unit_converter(breaks_x))) + 
          geom_rect(aes(xmin = START_g, xmax = END_g, ymin = -2, ymax = -1.4), 
            data = g_region_dt, fill = 'red', colour = 'red', alpha = 0.5, size = 0.5, show.legend = FALSE)
    } else {
        max_xend = contig_df[nrow(contig_df), 'ends']
        p = p + 
          geom_rect(aes(xmin = GPOS_X, xmax = GPOS2_X, ymin = -2, ymax = -1.4, fill = stain), 
            data = b_gdf, colour = 'black', size = 0.1, show.legend = FALSE) + 
          geom_segment(aes(x = ends, y = -2.4, xend = ends, yend = -1), 
            data = head(contig_df, -1), color = 'gray80', show.legend = FALSE) + 
          coord_cartesian(xlim = c(0, max_xend), 
            ylim = c(-2.4, -1), expand = FALSE, clip = 'off') +
          scale_x_continuous(limits = c(0, max_xend + max_xend * .01), 
            breaks = contig_df[['centers']], labels = contig_df[['names']])
    }
    return(p)
}

### plot x-axis label
label_bottom = function(xlab) {
    p = ggplot() + 
      theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'), plot.background = element_blank(), 
        axis.title.x = element_text(size = 18, face = 'bold', colour = 'black', vjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    p = p + 
      labs(x = xlab, y = '', colour = 'black')
    return(p)
}

### generate layout of CN plots
#the distance edge and panel is 2 cm
CN_scatter = function(cn_gdf, cn_seg_gdf, sv_gdf, contig_df, region_df, cytoband_gdf, height_adj, xlab, ylab, dots_size) {
    # get the ylim
    cn_max = 2^ceiling(log2(max(cn_gdf[['CN']])))
    cn_height = (2 + (ceiling(log2(cn_max))-1)/2*2)*height_adj # the diff between log2(2^1) and log2(2^n) is (n-1)

    # plot CN in a scatter plot
    ## plot panel
    cnl_p = CN_point(cn_gdf=cn_gdf, cn_seg_gdf=cn_seg_gdf, contig_df=contig_df, cytoband_gdf=cytoband_gdf, cn_max=cn_max, dots_size=dots_size)
    gc <- ggplotGrob(cnl_p)
    gc$heights[unique(gc$layout[gc$layout$name == 'panel', 't'])] = unit(cn_height, 'cm')
    ## plot left y-axis (total margin is 1 + 1.5 + 1 = 3.5 cm)
    yl = yaxis_left(cn_gdf, cn_max, ylab)
    gyl <- ggplotGrob(yl)
    gyl$widths[unique(gyl$layout[gyl$layout$name == 'ylab-l', 'r'])] = unit(0.65, 'cm')
    gyl$widths[unique(gyl$layout[gyl$layout$name == 'axis-l', 'r'])] = unit(0.8, 'cm')
    gyl$widths[unique(gyl$layout[gyl$layout$name == 'panel', 'l'])] = unit(0.05, 'cm')
    gyl$heights[unique(gyl$layout[gyl$layout$name == 'panel', 't'])] = unit(cn_height, 'cm')
    ## Merge Branch into Master
    gc <- cbind(gyl, gc, size='first')
    gb=gc

    # plot SVs using bezier curves
    if (nrow(sv_gdf) > 0) {
        # plot CN in a scatter plot
        ## plot panel
        sva_p = SV_bezier(sv_gdf=sv_gdf, contig_df=contig_df, cn_max=cn_max)
        gs <- ggplotGrob(sva_p)
        gs$heights[unique(gs$layout[gs$layout$name == 'panel', 't'])] = unit(6, 'cm')
        ## plot left y-axis
        ys = yaxis0_left(cn_max, cn_max + 6)
        gys <- ggplotGrob(ys)
        ## Merge Branch into Master
        gs <- cbind(gys, gs, size='last')
        gb <- rbind(gs, gb, size='last')
    }

    # x-axis
    if ((nrow(contig_df) == 1)&(ncol(region_df) == 3)) {
        # plot the bottom x-axis
        ## plot x-axis
        xa = xaxis_bottom(contig_df=contig_df, region_df=region_df)
        gx <- ggplotGrob(xa)
        gx$heights[unique(gx$layout[gx$layout$name == 'panel', 't'])] = unit(0, 'cm')
        gx$heights[unique(gx$layout[gx$layout$name == 'axis-b', 't'])] = unit(0.5, 'cm')
        ## plot left y-axis
        yx = yaxis0_left(-1.5, -1)
        gyx <- ggplotGrob(yx)
        ## Merge Branch into Master
        gx <- cbind(gyx, gx, size='last')
        gb <- rbind(gb, gx, size='first')
    }

    # cytoband
    if (nrow(cytoband_gdf) > 0) {
        # plot cytoband in a rectangle
        ## plot panel
        cbr_p = cytoband_rect(cytoband_gdf=cytoband_gdf, contig_df=contig_df, g_region_dt=region_df)
        gr <- ggplotGrob(cbr_p)
        gr$heights[unique(gr$layout[gr$layout$name == 'panel', 't'])] = unit(1, 'cm')
        gr$heights[unique(gr$layout[gr$layout$name == 'axis-b', 't'])] = unit(0.8, 'cm')
        ## plot left y-axis
        yr = yaxis0_left(-2.8, -1)
        gyr <- ggplotGrob(yr)
        ## Merge Branch into Master
        gr <- cbind(gyr, gr, size='last')
        gb <- rbind(gb, gr, size='first')
    }

    # plot the bottom label
    ## plot label
    lb = label_bottom(xlab=xlab)
    gl <- ggplotGrob(lb)
    gl$heights[unique(gl$layout[gl$layout$name == 'panel', 't'])] = unit(0, 'cm')
    gl$heights[unique(gl$layout[gl$layout$name == 'xlab-b', 't'])] = unit(1, 'cm')
    ## plot left y-axis
    yb = yaxis0_left(-3.5, -3)
    gyb <- ggplotGrob(yb)
    ## Merge Branch into Master
    gl <- cbind(gyb, gl, size='last')
    gb <- rbind(gb, gl, size='first')

    # return the ggplotGrob
    return(gb)
}


## segmentation
### Vector slider http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
Vector_slider = function(in_vt, window_n, step_n=1) {
    in_len = length(in_vt)
    vs_start = seq(from = 1, to = (in_len - window_n + 1), by = step_n)
    out_vt = vector(length = length(vs_start))
    for(i in 1:length(vs_start)){
        out_vt[i] = median(in_vt[vs_start[i]:(vs_start[i] + window_n - 1)])
    }
    return(out_vt)
}

### CN stabilizer
CN_polisher = function(cn_vt, max_iter_num, bin_n, label_vt='CN') {
    # get input
    cn0_vt = cn_vt
    # get index
    slider_ix = seq(from = 1, to = (length(cn0_vt) - bin_n + 1), by = 1) + floor((bin_n)/2)
    # loop until the pattern is stable to remove some tiny events
    for (i in 1:max_iter_num) {
        cn1_vt = Vector_slider(in_vt=cn0_vt, window_n=bin_n)
        cn1_vt = round(cn1_vt)
        cnx_vt = cn0_vt
        cnx_vt[slider_ix] = cn1_vt
        cnx_vt[1:(floor((bin_n)/2))] = cn1_vt[1]
        cnx_vt[(length(cn0_vt) - floor((bin_n)/2) + 1):length(cn0_vt)] = cn1_vt[length(cn1_vt)]
        cn1_vt = cnx_vt
        if (isTRUE(all.equal(cn0_vt, cn1_vt))) {
            cat(paste0('CN pattern is stable in ', label_vt, '.\n'));
            break
        }
        cn0_vt = cn1_vt
    }
    return(cn1_vt)
}

### CN binder
CN_binder = function(scn_dt) {
    # determine change points
    scn_dt[, smooth0_cp := c(1, diff(smooth_cn))]
    scn_dt[, smooth1_cp := c(diff(smooth_cn), 1)]
    # call start and end
    scns_dt = scn_dt[(smooth0_cp != 0), .(CONTIG, START, smooth_cn)]
    scne_dt = scn_dt[(smooth1_cp != 0), .(END, smooth_cn)]
    colnames(scns_dt)[3] = 'start_cn'
    colnames(scne_dt)[2] = 'end_cn'
    scn_sdt = cbind(scns_dt, scne_dt)
    # validate output segmentation data
    scn_sdt[, diff_cn := abs(end_cn - start_cn)]
    if (sum(scn_sdt[, diff_cn]) > 0) stop("In segmentation table, found a segment with unequal CN at start and end")
    # output
    #scn_sdt[, CN := start_cn]
    scn_sdt = scn_sdt[, .(CONTIG, START, END)]
    return(scn_sdt)
}

CN_doseg = function(cnj_dt, bgj_dt, cn_cname='PH1_CN', bin_n=9, max_iter_num=10) {
    # get bin length
    bin_len = cnj_dt[1, END - START + 1]
    # segmentation
    colnames(cnj_dt)[grep(cn_cname, colnames(cnj_dt))[1]] = 'CN'
    cnj_cn = cnj_dt[, CN]
    round_cn = round(cnj_cn)
    region_label = paste0(bgj_dt[,chr], ':', bgj_dt[, start], '-', bgj_dt[, end], ':', gsub('_CN', '', cn_cname))
    smooth_cn = CN_polisher(cn_vt=round_cn, bin_n=bin_n, max_iter_num=max_iter_num, label_vt=region_label)
    scn_dt = cbind(cnj_dt, smooth_cn=smooth_cn)
    scn_sdt = CN_binder(scn_dt=scn_dt)
    setkey(scn_sdt, CONTIG, START, END)
    # recalculate average CN
    scn_sdt = foverlaps(scn_sdt, cnj_dt)
    scn_sdt = scn_sdt[, .(median(CN)), by = .(CONTIG, i.START, i.END)]
    # build a new CN table
    colnames(scn_sdt) = c('CONTIG',  'START', 'END', 'CN')
    scn_sdt[, BINS := (END - START + 1) / bin_len]
    scn_sdt = cbind(data.table(SAMPLE=gsub('_CN', '', cn_cname)), scn_sdt)
    setcolorder(scn_sdt, c('SAMPLE', 'CONTIG',  'START', 'END', 'BINS', 'CN'))
    # output
    return(scn_sdt)
}

### Flipped CN segmentation
CN_segmentation = function(cn_dt, b_dt, bin_n=9, max_iter_num=10) {
    cn_sdt = data.table()
    # annotate G-banding regions
    bg_dt = copy(b_dt)
    bg_dt[, is_gband := ifelse((stain %like% '^g'), 1, 0)]
    # segmentation
    for (i in unique(cn_dt[, CONTIG])) {
        # get tables for current chr
        cni_dt = cn_dt[CONTIG == i]
        bgi_dt = bg_dt[chr == i]
        # split cytoBand table by non G-banding regions
        bgi_dt[, gband0_cp := c(1, diff(is_gband))]
        bgi_dt[, gband1_cp := c(diff(is_gband), 1)]
        bgis_dt = bgi_dt[(is_gband == 1) & (gband0_cp != 0), .(chr, start, stain)]
        bgie_dt = bgi_dt[(is_gband == 1) & (gband1_cp != 0), .(end, stain)]
        # get filtered table for G-banding regions
        colnames(bgis_dt)[3] = 'start_stain'
        colnames(bgie_dt)[2] = 'end_stain'
        bgif_dt = cbind(bgis_dt, bgie_dt)
        setkey(bgif_dt, chr, start, end)
        # validate filtered table for G-banding regions
        bgif_dt[, diff_stain := ifelse(start_stain != start_stain, 1, 0)]
        if (sum(bgif_dt[, diff_stain]) > 0) stop("In cytoband table, found a segment with different stain at start and end")
        # segmentation for each G-banding region
        for (j in 1:nrow(bgif_dt)) {
            bgj_dt = bgif_dt[j]
            cnj_dt = foverlaps(bgj_dt, cni_dt)
            colnames(cnj_dt)[grep('chr', colnames(cnj_dt))[1]] = 'CONTIG'
            setkey(cnj_dt, CONTIG, START, END)
            if (nrow(foverlaps(bgj_dt, cni_dt, nomatch=NULL)) < bin_n) next
            # do segmentation for H1
            h1_sdt = CN_doseg(cnj_dt=cnj_dt, bgj_dt=bgj_dt, cn_cname='PH1_CN', bin_n=bin_n, max_iter_num=max_iter_num)
            cn_sdt = rbind(cn_sdt, h1_sdt)
            if ('PH2_CN' %in% colnames(cn_dt)) {
              h2_sdt = CN_doseg(cnj_dt=cnj_dt, bgj_dt=bgj_dt, cn_cname='PH2_CN', bin_n=bin_n, max_iter_num=max_iter_num)
              cn_sdt = rbind(cn_sdt, h2_sdt)
            }
        }
    }
    setkey(cn_sdt, SAMPLE, CONTIG, START, END)
    return(cn_sdt)
}

### generate contour plot
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
        axis.title = element_text(size = 18, face = 'bold'), axis.text = element_text(size = 18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
      labs(x = gsub('_', ' ', cn4), y = gsub('_', ' ', cn5))
    return(p)
}

### generate layout of contour plots
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
     help = "Structural variation data from plotACN. (optinal) [%default]", metavar = "character"),
    make_option(c("-b", "--cytoband_tsv"), type = "character", default = "",
     help = "Cytoband data in TSV format. Please refer to 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz' (optinal) [%default]", metavar = "character"),
    make_option(c("-g", "--region"), type = "character", default = "",
     help = "The output region [%default] \n\t\te.g. g='chr17:36000000-42000000'", metavar = "character"),
    make_option(c("-m", "--cn_seg"), type = "character", default = "auto",
     help = "Whether perform segmentation or provide a seg file [%default]", metavar = "character"),
    make_option(c("--xlab"), type = "character", default = "",
     help = "x-axis label for output [%default]", metavar = "character"),
    make_option(c("--ylab"), type = "character", default = "Allelic CN",
     help = "y-axis label for output [%default]", metavar = "character"),
    make_option(c("--th_nite"), type = "double", default = "10",
     help = "Threshold for number of iterations [%default]", metavar = "double"),
    make_option(c("--th_nbin"), type = "double", default = "9",
     help = "Threshold for number of bins (roughly, 3 for 5kb, 5 for 10kb, 9 for 25kb) [%default]", metavar = "double"),
    make_option(c("--plot_width"), type = "double", default = "12",
     help = "Width of the ouput plot [%default]", metavar = "double"),
    make_option(c("--plot_height"), type = "double", default = "7",
     help = "Height of the ouput plot [%default]", metavar = "double"),
    make_option(c("--height_adj"), type = "double", default = "1",
     help = "Factor to adjust the height of each CN level [%default]", metavar = "double"),
    make_option(c("--dots_color"), type = "character", default = "black,gray",
     help = "Color of the dots in scatter plot [%default] \n\t\tIt can be one color for both 'black', or two color for each 'black,gray'", metavar = "character"),
    make_option(c("--dots_shape"), type = "character", default = "circle open,circle filled",
     help = "Shape of the dots in scatter plot [%default]", metavar = "character"),
    make_option(c("--dots_size"), type = "double", default = "1",
     help = "Size of the dots in scatter plot [%default]", metavar = "double"),
    make_option(c("--segs_color"), type = "character", default = "red,blue",
     help = "Color of the segments in scatter plot [%default] \n\t\tIt can be one color for both 'black', or two color for each 'black,gray'", metavar = "character"),
    make_option(c("--contour"), type="logical", default="FALSE", action = "store_true",
     help="Should a contour plot be generated [%default].", metavar="logical"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "plotFCN",
     help = "Output prefix [%default]", metavar = "character"))

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_number.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The following columns must be included: CONTIG, START, END and LOG2_COPY_RATIO.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

s = opt$options$sample_name
d = opt$options$ref_dict
k = opt$options$gene_bed
n = opt$options$sv_tsv
b = opt$options$cytoband_tsv
g = opt$options$region
x = opt$options$neu_cn
o = opt$options$output_prefix
m = opt$options$cn_seg
xl = opt$options$xlab
yl = opt$options$ylab
thi = opt$options$th_nite
thb = opt$options$th_nbin
pw = opt$options$plot_width
ph = opt$options$plot_height
chj = opt$options$height_adj
doc = opt$options$dots_color
dos = opt$options$dots_shape
doz = opt$options$dots_size
sgc = opt$options$segs_color
ct = opt$options$contour

f = opt$args[1]



# main
## read
### path
rv0 = as.numeric(R.version$major)
rv1 = as.numeric(R.version$minor)
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    f = ifelse(f == '-', 'file:///dev/stdin', f)
    n = ifelse(n == '-', 'file:///dev/stdin', n)
    if (file.exists(b)) {if(grepl('\\.gz$', b)) b_dt = fread(cmd = paste0('zcat ', b), skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
} else {
    f = ifelse(f == '-', 'cat /dev/stdin', f)
    n = ifelse(n == '-', 'cat /dev/stdin', n)
    if (file.exists(b)) {if(grepl('\\.gz$', b)) b_dt = fread(paste0('zcat ', b), skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
}
### fread
cn_dt = fread(f, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(n)) sv_dt = fread(n, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(b)) {if(!grepl('\\.gz$', b)) b_dt = fread(b, skip = 'chr', sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)}
if (file.exists(d)) contig_dict_df = read.delim(d, header=FALSE)
if (file_non0(k)) gene_dt = fread(k, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = FALSE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file_non0(m)) seg_dt = fread(m, sep = "\t", colClasses=list(character=1), stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)


## format
### process Copy Number data
cn_dt = CN_extractor(cn_dt=cn_dt)
### process Structural Variation data
if (file.exists(n)) sv_dt = SV_extractor(sv_dt=sv_dt)
### process the columns in cytoBand table
if (file.exists(b)) b_dt = CB_extractor(b_dt=b_dt)
### clip Genomic Region
if (grepl('[0-9XY]+', g)) {
    g_dt = GR_extractor(g_region=g)
    cn_dt = GR_clipper(data_dt=cn_dt, ref_dt=g_dt)
} else {
    g_dt = data.table()
}
### clip Gene Table
if (exists('gene_dt')) {
    gene_dt = GT_extractor(data_dt=cn_dt, gene_dt=gene_dt)
} else {
    gene_dt = data.table()
}
### clip Gene Table
if (exists('seg_dt')) {
    seg_dt = SG_extractor(seg_dt=seg_dt)
    if(grepl('[0-9XY]+', g)) seg_dt = SG_clipper(data_dt=seg_dt, ref_dt=g_dt)
} else {
    seg_dt = data.table()
}
### format dots color
doc_vt = unlist(strsplit(doc, ','))
map['PH1_color'] = doc_vt[1]
map['PH2_color'] = ifelse(length(doc_vt) == 1, doc_vt[1], doc_vt[2])
### format dots shape
dos_vt = unlist(strsplit(dos, ','))
map['PH1_shape'] = dos_vt[1]
map['PH2_shape'] = ifelse(length(dos_vt) == 1, dos_vt[1], dos_vt[2])
### format segs color
sgc_vt = unlist(strsplit(sgc, ','))
map['PH1'] = sgc_vt[1]
map['PH2'] = ifelse(length(sgc_vt) == 1, sgc_vt[1], sgc_vt[2])


## calculate
### CN segmentation
if (nrow(seg_dt) > 0) {
    cn_seg_dt = seg_dt
} else if (m == 'auto') {
    cn_seg_dt = CN_segmentation(cn_dt=cn_dt, b_dt=b_dt, bin_n=thb, max_iter_num=thi)
} else {
    cn_seg_dt = data.table()
}
### CN calling
cn_dt = CN_caller(cn_dt, gene_dt)

## map
### assemble contig
if (!file.exists(d)) contig_dict_df = data.table()
contig_dt = contig_assembler(ref_dict_df=contig_dict_df, g_region_dt=g_dt, cn_dt=cn_dt)
### map cn_dt to contig
cn_gdt = cn2contig(cn_dt=cn_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### map cn_sdt to contig
cn_seg_gdt = data.table()
if (nrow(cn_seg_dt) > 0) cn_seg_gdt = cns2contig(cn_seg=cn_seg_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### map sv_dt to contig
if (!file.exists(n)) sv_gdt = data.table()
if (file.exists(n)) sv_gdt = sv2contig(sv_dt=sv_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### map cytoband to contig
if (!file.exists(b)) b_gdt = data.table()
if (file.exists(b)) b_gdt = cb2contig(b_dt=b_dt, contig_dt=contig_dt, g_region_dt=g_dt)
### convert to data.frame
g_df = as.data.frame(g_dt)
if (file.exists(n)) sv_df = as.data.frame(sv_dt)
contig_df = as.data.frame(contig_dt)
cn_gdf = as.data.frame(cn_gdt)
cn_seg_gdf = as.data.frame(cn_seg_gdt)
sv_gdf = as.data.frame(sv_gdt)
b_gdf = as.data.frame(b_gdt)


## plot
### print scatter plot for CN
pdf(paste0(o, '.scatter.pdf'), pw, ph, useDingbats=FALSE)
    cs_gb=CN_scatter(cn_gdf=cn_gdf, cn_seg_gdf=cn_seg_gdf, sv_gdf=sv_gdf, contig_df=contig_df, region_df=g_df, cytoband_gdf=b_gdf, height_adj=chj, xlab=xl, ylab=yl, dots_size=doz)
    grid.draw(cs_gb)
dev.off()
### print contour plot for CN
if (ct & ('PH2_CN' %in% colnames(cn_dt))) {
    cn_df = as.data.frame(cn_dt)
    cn4_max = 2^ceiling(log2(max(cn_df[, 4])))
    cn5_max = 2^ceiling(log2(max(cn_df[, 5])))
    pdf(paste0(o, '.contour.pdf'), useDingbats=FALSE, 
      width = (3 + (ceiling(log2(cn4_max))-1)/2*2) / 2.54 * 1.5, 
      height = (3 + (ceiling(log2(cn5_max))-1)/2*2) / 2.54 * 1.5)
        cc_gb=CN_scatter2(cn_df=cn_df, xmax=cn4_max, ymax=cn5_max)
        grid.draw(cc_gb)
    dev.off()
}


## write
### print cn_gdf
colnames(cn_gdf)[grep('CN_type', colnames(cn_gdf))] = 'SAMPLE'
cn_gdf[['SAMPLE']] = paste0(s, '.', cn_gdf[['SAMPLE']])
cn_gdf = cn_gdf[, c('SAMPLE', 'CONTIG', 'START', 'END', 'CN')]
cn_gdf[['SAMPLE']] = gsub('black', 'PH1.black', cn_gdf[['SAMPLE']])
cn_gdf[['SAMPLE']] = gsub('gray', 'PH2.gray', cn_gdf[['SAMPLE']])
write.table(cn_gdf, paste0(o, '.CN.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
### print cn_seg_gdf
if (nrow(cn_seg_gdf) > 0) {
    cn_seg_gdf = cn_seg_gdf[, c('SAMPLE', 'CONTIG', 'START', 'END', 'BINS', 'CN')]
    cn_seg_gdf[['SAMPLE']] = paste0(s, '.', cn_seg_gdf[['SAMPLE']])
    cn_seg_gdf[['SAMPLE']] = gsub('red', 'PH1.red', cn_seg_gdf[['SAMPLE']])
    cn_seg_gdf[['SAMPLE']] = gsub('blue', 'PH2.blue', cn_seg_gdf[['SAMPLE']])
    write.table(cn_seg_gdf, paste0(o, '.CN.seg'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
}


# test
#s = 'SM_9_11_EAC1'
#f = '/home/unix/cbao/tmp/cn_phased_9_11-EAC1_chr17_nov14.dat'
#n = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/SV_evolution/CZ_bkps/EAC_9_11_bkps.txt'
#d = '/cga/bass/Chunyang/ref/hg19/Homo_sapiens_assembly19.1_22.dict'
#k = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/PhylogicNDT/Mutect2_gga_filtered1GF_forecalled_census_bed/EAC-9_11-EAC1.census.mut_cnv.bed'
#p = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p1p2/het_allelic_counts_normal_t0_phased1f_vcf/EAC-9_11-NORM.hets.phased1f.vcf.gz'
#b = '/cga/bass/Chunyang/ref/hg19/cytoBand.txt'
#o = '/home/unix/cbao/tmp/EAC-9_11-EAC1.chr17.5kb'
#g = 'chr17:35000000-45000000'
