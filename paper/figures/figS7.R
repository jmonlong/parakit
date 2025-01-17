suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(rjson))

## function to split reads/paths when they loop back in the pangenome
splitPaths <- function(df){
  ## if read looks like it's aligning in reverse, flip
  ## message(df$read[1])
  if(mean(diff(df$node)<0)>.5){
    if(any(colnames(df) == 'readpos')){
      df = arrange(df, desc(readpos))
    }
    if(any(colnames(df) == 'ppos')){
      df = arrange(df, desc(ppos))
    }
  }
  ## assumes the input is sorted by position in the path/read
  ## is passing through the cycling node
  cyc.pos = which(df$class=='cyc_r')
  ## keep only if the cycling is not near the end of the read
  cyc.pos = cyc.pos[which(cyc.pos/nrow(df)>.1 & cyc.pos/nrow(df)<.9)]
  ## if no need to split, return same input
  if(length(cyc.pos)==0) {
    df$path_part = 1
    return(df)
  }
  ## otherwise split each interval
  cyc.pos = c(1, cyc.pos, nrow(df))
  df = lapply(1:(length(cyc.pos)-1), function(ii){
    dfs = df[cyc.pos[ii]:cyc.pos[ii+1],]
    dfs$path_part = ii
    dfs
  }) %>% bind_rows
  return(df)
}

## which "scale" to use, either 'pangenome' (default) or 'genome'
fig.scale = 'pangenome'
## how much to shift the nodes from module 1/2 down/up
## (between 0 and 0.5 to not overlap with other parts of the graph)
path.v.shift.reads = .2
## how much to vertically shift the points of different classes (for aesthetic purpose)
path.v.shift.haps = 20

samples = c('DEN63190')
## genes.to.show = c('CYP21A1P', 'CYP21A2')
genes.to.show = c("C4A", "C4B", "TNXA", "TNXB", "CYP21A2", "CYP21A1P")

## offset for this region
config = fromJSON(file='../data/rccx.grch38_hprc.mc.config.json')
c1_s = as.numeric(gsub('.+:(.+)-.+', '\\1', config$c1))
c1_e = as.numeric(gsub('.+:.+-(.+)', '\\1', config$c1))
reg.offset = 1 + c1_s - config$flank_size

## load node info
ninfo = read.table('../data/rccx.grch38_hprc.mc.node_info.tsv', as.is=TRUE, header=TRUE)
ninfo = ninfo %>% filter(ref + c1 + c2 > 0) %>% arrange(node)

## assign a position on the second module for each node
ninfo$pos = NA
last_pos = 0
for(ii in 1:nrow(ninfo)){
  if(ninfo$rpos_min[ii] == ninfo$rpos_max[ii]){
    ninfo$pos[ii] = last_pos
  } else {
    ninfo$pos[ii] = ninfo$rpos_max[ii]
    last_pos = ninfo$rpos_max[ii]
  }
}
ninfo$pos[which(ninfo$pos==0)] = sort(unique(ninfo$pos))[2]
ninfo$pos = ninfo$pos + reg.offset

## load gene elements annotation
g.df = read.table("../data/CYP21A2.gencodev43.nearby_genes.tsv", as.is=TRUE, sep='\t', header=TRUE)
g.df = subset(g.df, gene_name %in% genes.to.show)

## load data
reads.df = lapply(samples, function(samp){
  read.table(paste0('../data/', samp, '.all_small.pdf.tsv'), as.is=TRUE, header=TRUE, comment.char="") %>%
    mutate(sample=samp)
}) %>% bind_rows %>% merge(ninfo)

vars = lapply(samples, function(samp){
  read.table(paste0('../data/', samp, '.calls.tsv'), as.is=TRUE, header=TRUE, check.names=F) %>%
    mutate(sample=samp)
}) %>% bind_rows %>% mutate(variant=factor(variant, unique(variant)),
                            allele=factor(allele, c('alt', 'ref', 'NA')))

reads.i = reads.df %>% filter(read %in% vars$read)
## pick longest read per variant downsample the reads to show
rl.df = reads.i %>% group_by(read) %>% summarize(length=n())
reads.to.plot = vars %>% filter(allele=='alt') %>%
  select(sample, variant, read) %>% unique %>% 
  merge(rl.df) %>%
  arrange(desc(length)) %>% group_by(sample, variant) %>% do(head(., 1)) %>% .$read
reads.i = reads.i %>% filter(read %in% as.character(reads.to.plot))

## load haplotypes
haps = lapply(samples, function(samp){
  stats = read.table(paste0('../data/', samp, '.paths-stats.tsv'), as.is=TRUE, header=TRUE)
  haps.i = read.table(paste0('../data/', samp, '.paths-info.tsv'), as.is=TRUE, header=TRUE)
  ## if same haplotype selected twice, duplicate with a different name
  if(stats$hap1[1] == stats$hap2[1]){
    new_hapn = paste0(stats$hap1[1], '_d')
    stats$hap2[1] = new_hapn
    haps.i = haps.i %>% filter(hap==stats$hap1[1]) %>%
      mutate(hap=new_hapn) %>% rbind(haps.i)
  }
  ## keep info for the best pair only
  haps.i %>% filter(hap %in% c(stats$hap1[1], stats$hap2[1])) %>% mutate(sample=samp)
}) %>% bind_rows

## ## Genomic element annotation
## project on nodes within the collapsed part of the pangenome, i.e. between cycling nodes
cyc_nodes = sort(subset(ninfo, class %in% c('cyc_l', 'cyc_r'))$node)
ninfo.col = ninfo %>% filter(node >= cyc_nodes[1], node <= cyc_nodes[2])
## make GRange objects for nodes and gene elements
ninfo.min.gr = GRanges(g.df$chr[1], IRanges(ninfo.col$rpos_min, width=ninfo.col$size))
ninfo.max.gr = GRanges(g.df$chr[1], IRanges(ninfo.col$rpos_max, width=ninfo.col$size))
g.gr = g.df %>% mutate(start=start-reg.offset, end=end-reg.offset) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)
## match them
ol = rbind(
  findOverlaps(g.gr, ninfo.min.gr) %>% as.data.frame %>%
  mutate(node=ninfo.col$node[subjectHits]),
  findOverlaps(g.gr, ninfo.max.gr) %>% as.data.frame %>%
  mutate(node=ninfo.col$node[subjectHits])) %>% 
  group_by(queryHits) %>% summarize(nstart=min(node), nend=max(node))
## add node start/end information to the gene annotation
g.df$nstart = NA
g.df$nstart[ol$queryHits] = ol$nstart
g.df$nend = NA
g.df$nend[ol$queryHits] = ifelse(ol$nstart==ol$nend, ol$nend + 1, ol$nend)
## annotate genes as being in module 1 or 2
c1c2.lim = ninfo %>% filter(class=='cyc_l') %>% .$rpos_max
g.df = g.df %>% mutate(module=ifelse((start+end)/2-reg.offset<c1c2.lim, 'c1', 'c2'))
## if we're plotting in "genome" scale, we need to shift the genes in module 1, annoying...
## let's begin by shifting by the median shift across the region
shift.bp = ninfo %>% filter(rpos_min != rpos_max) %>% mutate(shift=rpos_max-rpos_min) %>% .$shift %>% median
shift.lim = min(ninfo$pos)
g.df = g.df %>% mutate(shift=(start+end)/2 < shift.lim,
                       pstart=ifelse(shift, start + shift.bp, start),
                       pend=ifelse(shift, end + shift.bp, end))
## prepare TSS info
tss.df = g.df %>% filter(type=='gene') %>% group_by(gene_name, module) %>%
  summarize(node=ifelse(strand=='+', nstart, nend),
            pos=ifelse(strand=='+', pstart, pend), .groups='drop')


## graph list to store the panels of the final graph
ggp = list()
vlines.pos = c()
for(samp in samples){
  ## ## prepare read graphs
  ggp.vars.df = reads.i %>%
    filter(sample==samp, read %in% as.character(reads.to.plot)) %>%
    group_by(read) %>% arrange(readpos) %>% do(splitPaths(.)) %>%
    merge(., subset(vars, sample==samp & allele=='alt'), all.x=TRUE) %>% 
    mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
           ntype=ifelse(class!='both', 'module-specific', 'shared')) %>% 
    filter(!is.na(class))
  ## the points to highlight the variant calls on the reads
  ggp.vars.pts = ggp.vars.df %>% 
    group_by(read) %>%
    filter(!is.na(variant) & path_part == max(path_part)) %>% 
    mutate(vartype=ifelse(!is.na(node_u), 'fus', 'snp'),
           vartype=factor(vartype, c('fus', 'snp')),
           node=ifelse(!is.na(node_u), (node + node_u)/2, node),
           variant=ifelse(!is.na(node_u),
                          paste0(round((pos_l_1 + pos_u_1)/2), '_FUS_',
                                 round((pos_l_2 + pos_u_2)/2)), as.character(variant)))
  vlines.pos = c(vlines.pos, ggp.vars.pts$node)
  ## rename reads
  reads.o = unique(ggp.vars.df$read)
  ggp.vars.pts = ggp.vars.pts %>%
    mutate(read=factor(read, levels=reads.o,
                       labels=paste0('supp. read ', 1:length(unique(reads.o)))))
  ggp.vars.df = ggp.vars.df %>%
    mutate(read=factor(read, levels=reads.o,
                       labels=paste0('supp. read ', 1:length(unique(reads.o)))))
  ## prepare ggplot object
  ggp.o = ggplot(ggp.vars.df, aes(x=node, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) + xlab('node position in the pangenome')
  if(fig.scale == 'genome') {
    ggp$reads = ggplot(ggp.vars.df, aes(x=pos, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) + xlab('position in the chromosome')
  }  
  ggp.o = ggp.o + 
    geom_point(aes(shape=vartype), data=ggp.vars.pts, size=3) +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_shape_manual(values=2:1, breaks=c('fus', 'snp'), name='variant type') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('read\npart') +
    theme_bw() + 
    facet_grid(read~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=8, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))
  ggp[[paste0('reads_', samp)]] = ggp.o
  ## ## prepare allelic balance
  ## count the coverage of each position in the reads
  ## keeping only sites where both modules are supported
  ## if no support for module 2 but partial support for module 1, use 1-mod1_support
  reads.counts = reads.df %>% arrange(node) %>%
    filter(class %in% c('c1', 'c2'), sample==samp) %>% 
    group_by(rnode) %>% 
    summarize(site=ifelse(any('c1' %in% class), 'c1', NA),
              site=ifelse(any('c2' %in% class), 'c2', site),
              site=ifelse(all(c('c1', 'c2') %in% class), 'c1c2', site),
              c2.prop=mean(class=='c2'),
              c1.prop=mean(class=='c1'),
              c2.prop.adj=ifelse(c2.prop==0, 1-c1.prop, c2.prop),
              pos=pos[1], node=node[1], depth=n()) %>%
    filter(!is.na(site), site=='c1c2')
  ## ggplot object
  if(fig.scale == 'pangenome'){
    ggp.o = ggplot(reads.counts, aes(x=node, y=c2.prop.adj)) +
    xlab('node in collapsed pangenome')
  } else {
    ggp.o = ggplot(reads.counts, aes(x=pos, y=c2.prop.adj)) +
    xlab('position in chromosome')
  }
  ggp.o = ggp.o +
    geom_hline(yintercept=c(.25, .5, .75), color='#720e21', linetype=1, alpha=.5) +
    geom_hline(yintercept=c(1/3, 2/3), color='#03195b', linetype=1, alpha=.5) +
    scale_y_continuous(breaks=seq(0, 1, .25),
                       minor_breaks=c(),
                       limits=c(0,1)) + 
    geom_point(alpha=.7) + theme_bw() +
    ylab('module 2\nallelic ratio') +
    facet_grid("allele\nsupport"~.) +
    theme(legend.position='top', strip.text.y=element_text(angle=0))
  ggp[[paste0('as_', samp)]] = ggp.o
  ## ## coverage on non-specific nodes
  ## tile the nodes in bins of 100 bp
  ## (some nodes, e.g. flanking nodes, are huge, so we limite the size bias by splitting them)
  n.gr = reads.df %>% filter(class %in% c('ref', 'none')) %>%
    mutate(seqnames=node) %>%
    group_by(seqnames) %>% summarize(start=min(startpos), end=max(endpos)) %>% 
    select(seqnames, start, end) %>% makeGRangesFromDataFrame
  t.gr = unlist(tile(n.gr, width=10))
  ## prepare the reads "ranges" to compare with the tiled nodes
  r.gr = reads.df %>% filter(class %in% c('ref', 'none')) %>%
    mutate(seqnames=node, start=startpos, end=endpos) %>%
    select(seqnames, start, end) %>% makeGRangesFromDataFrame
  ## count how many reads overlap each node tile
  t.gr$cov = countOverlaps(t.gr, r.gr)
  ## convert tiles to data.frame and add node information
  cov.df = as.data.frame(t.gr) %>% mutate(node=as.integer(as.character(seqnames))) %>%
    select(-strand, -seqnames)
  cov.df$bin = factor(as.character(t.gr), levels=as.character(t.gr))
  cov.df = merge(cov.df, ninfo)
  ## assumes diploidy and normalize with coverage on the flanking region
  cov.df = cov.df %>% ungroup %>% mutate(cov.n=2*cov/median(cov[which(class=='ref')]))
  ## normalize using the edges flanking the module
  countFlankEdges <- function(read.df, ref.nodes, window.size=20){
    res_l = sapply(which(read.df$class=='cyc_l'), function(pos){
      up.pos = max(1, pos - window.size)
      if(any(ref.nodes %in% read.df$node[up.pos:pos])){
        return('fl')
      } else {
        return('cyc')
      }
    })
    res_r = sapply(which(read.df$class=='cyc_r'), function(pos){
      dw.pos = min(nrow(read.df), pos + window.size)
      if(any(ref.nodes %in% read.df$node[pos:dw.pos])){
        return('fl')
      } else {
        return('cyc')
      }
    })
    res = c(unlist(res_l), unlist(res_r))
    if(length(res) == 0) return(tibble(edge=NA))
    return(tibble(edge=res))
  }  
  ref.nodes = reads.df %>% filter(rpos_min==rpos_max) %>% .$node %>% unique
  edge.counts = reads.df %>% arrange(read, startpos) %>%
    group_by(read) %>% do(countFlankEdges(., ref.nodes))
  fl_c = sum(edge.counts$edge == 'fl', na.rm=TRUE) / 2
  cyc_c = sum(edge.counts$edge == 'cyc', na.rm=TRUE) / 2
  fl_cyc_cn = 2 * (fl_c + cyc_c) / fl_c
  ## message('Estimated copy number: ', fl_cyc_cn)
  cov.df = cov.df %>% ungroup %>% mutate(cov.n.2=2*cov/fl_c)  
  ## summarize the coverage for each node (across all tiles),
  ## keeping only nodes in both modules (or very close)
  cov.s = cov.df %>%
    filter(c1 > c2 * .8, c2 > c1 * .8, class!='ref') %>% 
    group_by(pos, node, class) %>% summarize(cov.n=median(cov.n),
                                             cov.n.2=median(cov.n.2))
  ## median coverage to show on the graph
  med.cov.2 = cov.df %>%
    filter(c1 > c2 * .8, c2 > c1 * .8, class!='ref') %>% .$cov.n.2 %>% median
  if(fig.scale == 'pangenome'){
    ggp.o = ggplot(cov.s, aes(x=node, y=cov.n.2)) + 
      xlab('node in collapsed pangenome')
  } else {
    ggp.o = ggplot(cov.s, aes(x=pos, y=cov.n.2)) + 
      xlab('position in chromosome')
  }
  ggp.o = ggp.o +
    geom_point(alpha=.7) +
    theme_bw() +
    geom_hline(yintercept=4, linetype=1, alpha=.1) +
    geom_hline(yintercept=med.cov.2, linetype=2, alpha=.7) +
    geom_hline(yintercept=fl_cyc_cn, linetype=3, alpha=.7) +
    ylab('estimated\ncopy number') +
    facet_grid("coverage\n"~.) +
    scale_y_continuous(breaks=seq(0,10,2)) + 
    theme(legend.position='top', strip.text.y=element_text(angle=0))
  ggp[[paste0('cov_', samp)]] = ggp.o  
  ## ## haplotypes
  ggp.haps.s.df = haps %>% filter(sample==samp) %>% merge(ninfo) %>%
    group_by(hap) %>% arrange(ppos) %>% do(splitPaths(.)) %>%
    ungroup %>% 
    mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
           ntype=ifelse(class!='both', 'module-specific', 'shared'),
           hap=factor(hap, levels=unique(hap), labels=c('H1', 'H2'))) %>%
    filter(!is.na(ntype))
  if(fig.scale == 'pangenome'){
    ggp.o = ggplot(ggp.haps.s.df, aes(x=node,
                                      y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
      xlab('node in collapsed pangenome')
  } else {
    ggp.o = ggplot(ggp.haps.s.df, aes(x=pos,
                                           y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
      xlab('position in chromosome')
  }
  ggp.o = ggp.o +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('path\npart') + 
    theme_bw() + 
    facet_grid(hap~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=8, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))
  ggp[[paste0('hap_', samp)]] = ggp.o    
}

## start the ggplot object
ggp$genes = g.df %>% filter(type %in% c('exon', 'gene')) %>% 
  ggplot(aes(color=module))
## palette with just 1st/3rd values (module 1/2 in the full palette)
pal.set1 = brewer.pal(6, 'Set1')[c(1,3)]
## define the columns to use
if(fig.scale == 'pangenome'){
  ggp$genes = ggp$genes + 
    geom_segment(aes(x=nstart, xend=nend, y=0, yend=0, linewidth=type)) +
    xlab('node ID in the pangenome')
} else {
  ggp$genes = ggp$genes + 
    geom_segment(aes(x=pstart, xend=pend, y=0, yend=0, linewidth=type)) +
    xlab('position in the chromosome')
}
## add the rest of the ggplot elements
ggp$genes = ggp$genes + 
  facet_grid(gene_name~.) +
  scale_linewidth_manual(values=c(3, 1)) +
  scale_y_continuous(breaks=0:1) +
  scale_color_manual(values=pal.set1) + 
  theme_bw() +
  guides(color='none') +
  labs(linewidth=NULL) +
  ylab('gene\nannotation') + 
  theme(strip.text.y=element_text(angle=0),
        axis.text.y=element_blank(), 
        legend.position=c(.01,.01), legend.justification=c(0,0))
## add TSS symbol
if(fig.scale == 'pangenome'){
  ggp$genes = ggp$genes + geom_point(aes(x=node, y=0), size=2, data=tss.df,
                                     shape=4, alpha=.7)
} else {
  ggp$genes = ggp$genes + geom_point(aes(x=pos, y=0), size=2, data=tss.df,
                                     shape=4, alpha=.7)
}

## shortcut to disable x-axis legend (for all panels except the bottom one)
nox = theme(axis.text.x=element_blank(), axis.title.x=element_blank())
noleg = guides(alpha='none', color='none', shape='none', linewidth='none')
nobg = theme(plot.background=element_blank())

## vertical dotted lines
vlines = geom_vline(xintercept=unique(vlines.pos), linetype=3, alpha=.4)

## define the x-axis zoom
ggp.xlims = xlim(min(g.df$nstart, na.rm=TRUE) - 300, max(g.df$nend, na.rm=TRUE) + 300)
ggp.xlims = xlim(min(ninfo$node), max(ninfo$node))

names(ggp)

## put the panels together
pdf('figS7.raw.pdf', 7, 6)
plot_grid(ggp$reads_DEN63190 + ggp.xlims + nox + nobg + vlines,
          ggp$cov_DEN63190 + ggp.xlims + nox + nobg + noleg + vlines,
          ggp$as_DEN63190 + ggp.xlims + nox + nobg + noleg + vlines,
          ggp$hap_DEN63190 + ggp.xlims + nox + nobg + noleg + vlines,
          ggp$genes + ggp.xlims + noleg + nobg + vlines,
          ncol=1, align='v', rel_heights=c(3,1.5,1.5,1.5,3))
dev.off()
