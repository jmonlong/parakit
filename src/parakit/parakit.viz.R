suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(rjson))

## define arguments
args = list()
args$viz = list(arg='-v', desc='visualization mode')
args$nodes = list(arg='-n', desc='node information in TSV')
args$config = list(arg='-j', desc='config json')
args$reads = list(arg='-r', desc='input alignments in TSV from "parakit readcov"')
args$genes = list(arg='-e', desc='input genome element annotation TSV')
args$calls = list(arg='-c', desc='calls TSV')
args$nreads = list(arg='-m', val=3,
                   desc='maximum number of supporting reads to show in graph')
args$label = list(arg='-l', desc='a label to use as title of the graphs', val='')
args$hstats = list(arg='-d', desc='diplotype paths, sorted')
args$hpaths = list(arg='-p', desc='haplotype paths information')
args$out = list(arg='-o', desc='output PDF file', val='parakit.viz.pdf')

## parse arguments
args.i = commandArgs(TRUE)
arg.to.arg = names(args)
names(arg.to.arg) = as.character(sapply(args, function(l) l$arg))
ii = 1
show.usage = FALSE
while(ii <= length(args.i) & !show.usage){
  if(any(args.i[ii] == names(arg.to.arg))){
    args[[arg.to.arg[args.i[ii]]]]$val = args.i[ii+1]
    ii = ii + 1
  } else {
    if(args.i[ii] != '-h'){
      cat("\n\nError:", args.i[ii], ' argument unknown.\n', file = stderr())
    }
    show.usage = TRUE
  }
  ii = ii + 1
}
if(show.usage){
  cat("\nUsage:\n", file = stderr())
  lapply(args, function(ll) cat(paste(ll$arg, '\t', ll$desc, "\n"), file=stderr()))
  stop()
}

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

## graph list to store the panels of the final graph
ggp = list()
## to save the boundaries of the x-axis
xlims_v = c()
## if variants, will store the corresponding vertical lines to add to other graphs
var.vl = NULL

## load node info
ninfo = read.table(args$nodes$val, as.is=TRUE, header=TRUE)
ninfo = subset(ninfo, ref + c1 + c2 > 0)

## offset for this region
config = fromJSON(file=args$config$val)
c1_s = as.numeric(gsub('.+:(.+)-.+', '\\1', config$c1))
c1_e = as.numeric(gsub('.+:.+-(.+)', '\\1', config$c1))
reg.offset = 1 + c1_s - config$flank_size

## load read information
reads.df = NULL
if(args$viz$val %in% c('calls', 'all', 'all_small', 'allele_support', 'annotate')){
  reads.df = read.table(args$reads$val, as.is=TRUE, header=TRUE) %>%
    merge(ninfo)
}


##
## Calls and supporting reads
##

if(args$viz$val %in% c('calls', 'all', 'all_small')){
  ## load read-variants table
  vars = read.table(args$calls$val, as.is=TRUE, header=TRUE, check.names=F)
  vars = vars %>% mutate(variant=factor(variant, unique(variant)),
                         allele=factor(allele, c('alt', 'ref', 'na')))

  ## cluster reads
  vars.m = vars %>% mutate(allele=as.numeric(allele)) %>% select(-node) %>%
    unique %>% 
    pivot_wider(id_cols=read, names_from=variant, values_from=allele)
  hc.o = hclust(dist(vars.m[,-1]))
  read.ord = vars.m$read[hc.o$order]
  vars$read = factor(vars$read, read.ord)

  ## keep relevant information about reads of interest
  reads.i = reads.df %>% filter(read %in% vars$read)

  ## potentially downsample the reads to show
  rl.df = reads.i %>% group_by(read) %>% summarize(length=n())
  ## select the longest ones though
  reads.to.plot = vars %>% filter(allele=='alt') %>%
    select(variant, read) %>% unique %>% 
    merge(rl.df) %>%
    arrange(desc(length)) %>% group_by(variant) %>% do(head(., as.numeric(args$nreads$val))) %>% .$read
  
  ## vertical dotted line to help following the variants called
  var.vl = geom_vline(xintercept=unique(vars$node), linetype=3, linewidth=.3)

  ## prepare data.frame for graph
  ## keep reads to show, split, add variant info
  ggp.vars.df = reads.i %>%
    filter(read %in% as.character(reads.to.plot)) %>% 
    group_by(read) %>% arrange(readpos) %>% do(splitPaths(.)) %>%
    merge(., subset(vars, allele=='alt'), all.x=TRUE) %>% 
    mutate(read=factor(read, intersect(read.ord, reads.to.plot)),
           class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
           ntype=ifelse(class!='both', 'module-specific', 'shared')) %>% 
    filter(!is.na(class))

  ## the points to highlight the variant calls on the reads
  ggp.vars.pts = ggp.vars.df %>% 
    group_by(read) %>%
    filter(!is.na(variant) & path_part == max(path_part))

  ## how much to shift the nodes from module 1/2 down/up
  ## (between 0 and 0.5 to not overlap with other parts of the graph)
  path.v.shift.reads = .2
  ## prepare ggplot object
  ggp$reads = ggplot(ggp.vars.df, aes(x=node, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
    var.vl + 
    geom_point(aes(shape=variant), data=ggp.vars.pts, size=3) +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('read\npart') + xlab('node position in the pangenome') +
    theme_bw() + 
    facet_grid(read~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))

  ## save the x-axis boundaries for later
  xlims_v = c(xlims_v, ggp.vars.df$node)
}


##
## Genomic element annotation
##

## load gene elements annotation
g.df = read.table(args$genes$val, as.is=TRUE, sep='\t', header=TRUE)

## subset to specific genes if specified in the config
if('genes' %in% names(config)){
  g.df = subset(g.df, gene_name %in% config$genes)
}

## make GRange objects for nodes and gene elements
ninfo.min.gr = GRanges(g.df$chr[1], IRanges(ninfo$rpos_min, width=ninfo$size))
ninfo.max.gr = GRanges(g.df$chr[1], IRanges(ninfo$rpos_max, width=ninfo$size))
g.gr = g.df %>% mutate(start=start-reg.offset, end=end-reg.offset) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

## match them
ol = rbind(
  findOverlaps(g.gr, ninfo.min.gr) %>% as.data.frame %>%
  mutate(node=ninfo$node[subjectHits]),
  findOverlaps(g.gr, ninfo.max.gr) %>% as.data.frame %>%
  mutate(node=ninfo$node[subjectHits])) %>% 
  group_by(queryHits) %>% summarize(nstart=min(node), nend=max(node))

## add node start/end information to the gene annotation
g.df$nstart = NA
g.df$nstart[ol$queryHits] = ol$nstart
g.df$nend = NA
g.df$nend[ol$queryHits] = ol$nend

## annotate genes as being in module 1 or 2
c1c2.lim = ninfo %>% filter(class=='cyc_l') %>% .$rpos_max
g.df = g.df %>% mutate(module=ifelse(start-reg.offset<c1c2.lim, 'c1', 'c2'))

## start the ggplot object
ggp$genes = g.df %>% filter(type %in% c('exon', 'gene')) %>% 
  ggplot(aes(color=module))

## if variants are also shown, add the vertical lines to help follow their positions
if(!is.null(var.vl)){
  ggp$genes = ggp$genes + var.vl
}

## palette with just 1st/3rd values (module 1/2 in the full palette)
pal.set1 = brewer.pal(6, 'Set1')[c(1,3)]

## add the rest of the ggplot elements
ggp$genes = ggp$genes + 
  geom_segment(aes(x=nstart, xend=nend, y=0, yend=0, linewidth=type)) + 
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
        legend.position=c(.01,.01), legend.justification=c(0,0)) +
  xlab('node ID in the pangenome')
tss.df = g.df %>% filter(type=='gene') %>% group_by(gene_name, module) %>%
  summarize(node=ifelse(strand=='+', nstart, nend), .groups='drop')
ggp$genes = ggp$genes + geom_point(aes(x=node, y=0), size=3, data=tss.df, shape=23, alpha=.7)

## save the x-axis boundaries for later
xlims_v = c(xlims_v, g.df$nstart, g.df$nend)

##
## allelic balance
##

if(args$viz$val %in% c('allele_support', 'all', 'all_small')){

  ## count the coverage of each position in the reads
  ## keeping only sites where both modules are supported
  ## if no support for module 2 but partial support for module 1, use 1-mod1_support
  reads.counts = reads.df %>% arrange(node) %>%
    group_by(rpos_min) %>% 
    summarize(site=ifelse(any('c1' %in% class), 'c1', NA),
              site=ifelse(any('c2' %in% class), 'c2', site),
              site=ifelse(all(c('c1', 'c2') %in% class), 'c1c2', site),
              c2.prop=mean(class=='c2'),
              c1.prop=mean(class=='c1'),
              c2.prop.adj=ifelse(c2.prop==0, 1-c1.prop, c2.prop),
              node=node[1], depth=n()) %>%
    filter(!is.na(site), site=='c1c2')

  ## ggplot object
  ggp$allele = ggplot(reads.counts, aes(x=node, y=c2.prop.adj)) +
    geom_point(alpha=.7) + theme_bw() +
    geom_hline(yintercept=.5, linetype=2) +
    geom_hline(yintercept=c(1/3,2/3), linetype=3) +
    ylim(0,1) +
    ylab('module 2\nallelic ratio') +
    xlab('node in collapsed pangenome') +
    facet_grid("allele\nsupport"~.) +
    theme(legend.position='top', strip.text.y=element_text(angle=0))

  ##
  ## coverage on non-specific nodes

  ## tile the nodes in bins of 100 bp
  ## (some nodes, e.g. flanking nodes, are huge, so we limite the size bias by splitting them)
  n.gr = reads.df %>% filter(class %in% c('ref', 'none')) %>%
    mutate(seqnames=node) %>%
    group_by(seqnames) %>% summarize(start=min(startpos), end=max(endpos)) %>% 
    select(seqnames, start, end) %>% makeGRangesFromDataFrame
  t.gr = unlist(tile(n.gr, width=100))

  ## in case we want to use reads with a minimum length, I leave this here
  rl.df = reads.df %>% group_by(read) %>% summarize(kbp=sum(endpos-startpos)/1e3)
  reads.sel = rl.df %>% filter(kbp>0) %>% .$read

  ## prepare the reads "ranges" to compare with the tiled nodes
  r.gr = reads.df %>% filter(read %in% reads.sel, class %in% c('ref', 'none')) %>%
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
  ## median coverage to show on the graph
  med.cov = cov.df %>% filter(class!='ref') %>% .$cov.n %>% median

  ## summarize the coverage for each node (across all tiles),
  ## keeping only nodes in both modules (or very close)
  cov.s = cov.df %>%
    filter(c1 > c2 * .9, c2 > c1 * .9, class!='ref') %>% 
    group_by(node, class) %>% summarize(cov.n=median(cov.n))

  ## ggplot object
  ggp$coverage = ggplot(cov.s, aes(x=node, y=cov.n)) +
    geom_point(alpha=.7) +
    theme_bw() +
    geom_hline(yintercept=4, linetype=2) +
    geom_hline(yintercept=med.cov, linetype=3) +
    ylab('estimated\ncopy number') +
    xlab('node in collapsed pangenome') +
    facet_grid("coverage\n"~.) +
    scale_y_continuous(breaks=seq(0,10,2)) + 
    theme(legend.position='top', strip.text.y=element_text(angle=0))
}

##
## haplotypes
##

if(args$viz$val %in% c('paths', 'all', 'all_small')){
  ## load stats on pairs of predicted haplotypes
  stats = read.table(args$hstats$val, as.is=TRUE, header=TRUE)
  ## load information for each predicted haplotype
  haps = read.table(args$hpaths$val, as.is=TRUE, header=TRUE)
  ## keep info for the best pair only
  haps = haps %>% filter(hap %in% c(stats$hap1[1], stats$hap2[1]))

  ## prep predicted haplotype paths
  ggp.haps.df = haps %>% mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
                                ntype=ifelse(class!='both', 'module-specific', 'shared'),
                                hap=factor(hap, levels=unique(hap), labels=c('hap_1', 'hap_2'))) %>% 
    filter(!is.na(class)) %>% arrange(class=='1', class=='2')

  ## how much to vertically shift the points of different classes (for aesthetic purpose)
  path.v.shift.haps = 20
  ## "dot plot" visualization: x=pangenome position, y=haplotype position
  ggp$haps = ggplot(ggp.haps.df, aes(x=node, y=ppos + path.v.shift.haps*(as.numeric(class) - 1),
                                     color=class, alpha=ntype)) +
    geom_point() +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.8,.05), name='node') + 
    facet_grid(hap~., scales='free', space='free') + 
    theme_bw() +
    theme(strip.text.y=element_text(angle=0),
          legend.position='top',
          axis.text.x=element_blank(), axis.title.x=element_blank()) +
    ylab('position in\npredicted haplotype')

  ## split paths visualization (like for reads above)
  ggp.haps.s.df = haps %>% merge(ninfo) %>%
    group_by(hap) %>% arrange(ppos) %>% do(splitPaths(.)) %>%
    ungroup %>% 
    mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
           ntype=ifelse(class!='both', 'module-specific', 'shared'),
           hap=factor(hap, levels=unique(hap), labels=c('hap_1', 'hap_2')))
  path.v.shift.reads = .2
  ggp$haps.s = ggplot(ggp.haps.s.df, aes(x=node, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('path\npart') + xlab('node position in the pangenome') +
    theme_bw() + 
    facet_grid(hap~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))
}

##
## annotate sequence
##

if(args$viz$val == 'annotate'){
  ## split paths visualization (like for reads above)
  ggp.annot.df = reads.df %>%
    group_by(read) %>% arrange(readpos) %>% do(splitPaths(.)) %>%
    ungroup %>% 
    mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
           ntype=ifelse(class!='both', 'module-specific', 'shared')) %>%
    filter(!is.na(class))
  path.v.shift.reads = .2
  ggp$annot = ggplot(ggp.annot.df,
                     aes(x=node, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('path\npart') + xlab('node position in the pangenome') +
    theme_bw() + 
    facet_grid(read~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))

  path.v.shift.reads = max(ggp.annot.df$readpos) * .02
  ggp$annot = ggplot(ggp.annot.df,
                     aes(x=node,
                         y=readpos + path.v.shift.reads*(as.numeric(class) - 1),
                         color=class, alpha=ntype)) +
    geom_point() +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.8,.05), name='node') + 
    facet_grid(read~., scales='free', space='free') + 
    theme_bw() +
    ylab('position in\ninput sequence') +
    xlab('node position in the pangenome') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))

  ## save the x-axis boundaries for later
  xlims_v = c(xlims_v, ggp.annot.df$node)
}

## to make sure all panels have the same x-axis limits 
ggp.xlims = xlim(min(xlims_v), max(xlims_v))

## add the vertical lines highlighting the variants calls (if present)
if(!is.null(var.vl) & 'coverage' %in% names(ggp)){
  ggp$coverage = ggp$coverage + var.vl
}
if(!is.null(var.vl) & 'allele' %in% names(ggp)){
  ggp$allele = ggp$allele + var.vl
}
if(!is.null(var.vl) & 'haps' %in% names(ggp)){
  ggp$haps = ggp$haps + var.vl
}
if(!is.null(var.vl) & 'haps.s' %in% names(ggp)){
  ggp$haps.s = ggp$haps.s + var.vl
}
if(!is.null(var.vl) & 'annot.s' %in% names(ggp)){
  ggp$annot.s = ggp$annot.s + var.vl
}
if(!is.null(var.vl) & 'annot' %in% names(ggp)){
  ggp$annot = ggp$annot + var.vl
}

## shortcut to disable x-axis legend (for all panels except the bottom one)
nox = theme(axis.text.x=element_blank(), axis.title.x=element_blank())

## put the panels together
pdf.h = ifelse(args$viz$val == 'all', 9, 6)
pdf(args$out$val, 9, pdf.h)

if(args$viz$val == 'allele_support'){
  plot_grid(ggp$coverage +
            ggp.xlims + nox,
            ggp$allele + 
            ggp.xlims + nox,
            ggp$genes + ggp.xlims + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(1,1,1))
}

if(args$viz$val == 'calls'){
  plot_grid(ggp$reads + 
            ggp.xlims + nox,
            ggp$genes + ggp.xlims + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(4,1.8))
}

if(args$viz$val == 'paths'){
  plot_grid(ggp$haps +
            ggp.xlims + nox,
            ggp$genes + ggp.xlims + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(2,1))
}

if(args$viz$val == 'annotate'){
  plot_grid(ggp$annot +
            ggp.xlims + nox,
            ggp$genes + ggp.xlims + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(2,1))
}

if(args$viz$val == 'all'){
  plot_grid(ggp$reads + ggp.xlims + nox,
            ggp$coverage + ggp.xlims + nox,
            ggp$allele + ggp.xlims + nox,
            ggp$haps.s + guides(alpha=FALSE, color=FALSE) + ggp.xlims + nox,
            ggp$genes + ggp.xlims + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(4,1.5,1.5,2, 2))
}

if(args$viz$val == 'all_small'){
  plot_grid(ggp$reads + ggp.xlims + nox + theme(legend.title=element_text(size=8),
                                                legend.text=element_text(size=8)),
            ggp$coverage + ggp.xlims + nox,
            ggp$allele + ggp.xlims + nox,
            ggp$haps.s + guides(alpha=FALSE, color=FALSE) + ggp.xlims + nox,
            ggp$genes + ggp.xlims  + labs(caption=args$label$val),
            ncol=1, align='v', rel_heights=c(2,1,1,1.2,2))
}

dev.off()
