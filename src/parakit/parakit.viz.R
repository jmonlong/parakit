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
args$scale = list(arg='-s', val='pangenome',
                   desc='pangenome or genome scale')
args$label = list(arg='-l', desc='a label to use as title of the graphs', val='')
args$hstats = list(arg='-d', desc='diplotype paths, sorted')
args$hpaths = list(arg='-p', desc='haplotype paths information')
args$out = list(arg='-o', desc='output PDF file', val='parakit.viz.pdf')

## parse arguments
args.i = commandArgs(TRUE)
## args.i = unlist(strsplit('-j rccx.grch38_hprc.mc.config.json -v all_small -n rccx.grch38_hprc.mc.node_info.tsv -e CYP21A2.gencodev43.nearby_genes.tsv -o results/rccx.grch38_hprc.mc/DEN63190/DEN63190.all_small.pdf -m 1 -c results/rccx.grch38_hprc.mc/DEN63190/DEN63190.calls.tsv -r results/rccx.grch38_hprc.mc/DEN63190/DEN63190.all_small.pdf.tsv -d results/rccx.grch38_hprc.mc/DEN63190/DEN63190.paths-stats.tsv -p results/rccx.grch38_hprc.mc/DEN63190/DEN63190.paths-info.tsv -l DEN63190', ' '))
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

## which "scale" to use, either 'pangenome' (default) or 'genome'
fig.scale = args$scale$val

## graph list to store the panels of the final graph
ggp = list()
## to save the boundaries of the x-axis
xlims_v = c()
## if variants, will store the corresponding vertical lines to add to other graphs
var.vl = NULL

## offset for this region
config = fromJSON(file=args$config$val)
c1_s = as.numeric(gsub('.+:(.+)-.+', '\\1', config$c1))
c1_e = as.numeric(gsub('.+:.+-(.+)', '\\1', config$c1))
reg.offset = 1 + c1_s - config$flank_size

## load node info
ninfo = read.table(args$nodes$val, as.is=TRUE, header=TRUE)
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

## min.fc = 2
## ninfo = ninfo %>%
##   mutate(class=ifelse(class=='none' & c2 > min.fc * c1, 'c2', class),
##          class=ifelse(class=='none' & c1 > min.fc * c2, 'c1', class))

## load read information
reads.df = NULL
if(args$viz$val %in% c('calls', 'all', 'all_small', 'allele_support', 'annotate')){
  reads.df = read.table(args$reads$val, as.is=TRUE, header=TRUE, comment.char="") %>%
    merge(ninfo)
}


##
## Calls and supporting reads
##

if(args$viz$val %in% c('calls', 'all', 'all_small')){

  ## load read-variants table
  vars = read.table(args$calls$val, as.is=TRUE, header=TRUE, check.names=F)
  vars = vars %>% mutate(variant=factor(variant, unique(variant)),
                         allele=factor(allele, c('alt', 'ref', 'NA')))

  ## cluster reads
  vars.m = vars %>% mutate(allele=as.numeric(allele)) %>% select(read, variant, allele) %>%
    mutate(allele=ifelse(is.na(allele), '0', allele)) %>% 
    unique %>% 
    pivot_wider(id_cols=read, names_from=variant, values_from=allele)
  read.ord = vars.m$read
  if(ncol(vars.m)>2){
    hc.o = hclust(dist(vars.m[,-1]))
    read.ord = vars.m$read[hc.o$order]
  }
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
    filter(!is.na(variant) & path_part == max(path_part)) %>% 
    mutate(node=ifelse(!is.na(node_u), (node + node_u)/2, node),
           variant=ifelse(!is.na(node_u),
                          paste0(round((pos_l_1 + pos_u_1)/2), '_FUS_',
                                 round((pos_l_2 + pos_u_2)/2)), as.character(variant)))

  ## vertical dotted line to help following the variants called
  var.vl = geom_vline(xintercept=unique(ggp.vars.pts$node), linetype=3, linewidth=.3)
  if(fig.scale == 'genome') {
    var.vl = geom_vline(xintercept=unique(ggp.vars.pts$pos), linetype=3, linewidth=.3)
  }
  
  ## how much to shift the nodes from module 1/2 down/up
  ## (between 0 and 0.5 to not overlap with other parts of the graph)
  path.v.shift.reads = .2
  ## prepare ggplot object
  ggp$reads = ggplot(ggp.vars.df, aes(x=node, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) + xlab('node position in the pangenome')
  if(fig.scale == 'genome') {
    ggp$reads = ggplot(ggp.vars.df, aes(x=pos, y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) + xlab('position in the chromosome')
  }
  
  ggp$reads = ggp$reads + var.vl + 
    geom_point(aes(shape=variant), data=ggp.vars.pts, size=3) +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('read\npart') +
    theme_bw() + 
    facet_grid(read~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))

  ## save the x-axis boundaries for later
  if(fig.scale == 'genome') {
    xlims_v = c(xlims_v, ggp.vars.df$pos)
  } else {
    xlims_v = c(xlims_v, ggp.vars.df$node)
  }
}

##
## allelic balance
##

if(args$viz$val %in% c('allele_support', 'all', 'all_small')){

  ## count the coverage of each position in the reads
  ## keeping only sites where both modules are supported
  ## if no support for module 2 but partial support for module 1, use 1-mod1_support
  reads.counts = reads.df %>% arrange(node) %>%
    filter(class %in% c('c1', 'c2')) %>% 
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
    ggp$allele = ggplot(reads.counts, aes(x=node, y=c2.prop.adj)) +
    xlab('node in collapsed pangenome')
  } else {
    ggp$allele = ggplot(reads.counts, aes(x=pos, y=c2.prop.adj)) +
    xlab('position in chromosome')
  }
  
  ggp$allele = ggp$allele +
    geom_hline(yintercept=c(.25, .5, .75), color='#720e21', linetype=1, alpha=.5) +
    geom_hline(yintercept=c(1/3, 2/3), color='#03195b', linetype=1, alpha=.5) +
    scale_y_continuous(breaks=seq(0, 1, .25),
                       minor_breaks=c(),
                       limits=c(0,1)) + 
    geom_point(alpha=.7) + theme_bw() +
    ylab('module 2\nallelic ratio') +
    facet_grid("allele\nsupport"~.) +
    theme(legend.position='top', strip.text.y=element_text(angle=0))
  
  ## save the x-axis boundaries for later
  if(fig.scale == 'genome') {
    xlims_v = c(xlims_v, reads.counts$pos)
  } else {
    xlims_v = c(xlims_v, reads.counts$node)
  }
  
  ##
  ## coverage on non-specific nodes

  ## tile the nodes in bins of 100 bp
  ## (some nodes, e.g. flanking nodes, are huge, so we limite the size bias by splitting them)
  n.gr = reads.df %>% filter(class %in% c('ref', 'none')) %>%
    mutate(seqnames=node) %>%
    group_by(seqnames) %>% summarize(start=min(startpos), end=max(endpos)) %>% 
    select(seqnames, start, end) %>% makeGRangesFromDataFrame
  t.gr = unlist(tile(n.gr, width=10))

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
  ## med.cov = cov.df %>% filter(class!='ref') %>% .$cov.n %>% median

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
  med.cov = cov.df %>%
    filter(c1 > c2 * .8, c2 > c1 * .8, class!='ref') %>% .$cov.n %>% median
  med.cov.2 = cov.df %>%
    filter(c1 > c2 * .8, c2 > c1 * .8, class!='ref') %>% .$cov.n.2 %>% median
  
  ## ggplot object
  ## ggp$coverage = ggplot(cov.s, aes(x=node, y=cov.n)) +
  ##   geom_point(alpha=.7) +
  ##   theme_bw() +
  ##   geom_hline(yintercept=4, linetype=2) +
  ##   geom_hline(yintercept=med.cov, linetype=3) +
  ##   geom_hline(yintercept=med.cov.2, linetype=4) +
  ##   geom_hline(yintercept=fl_cyc_cn, linetype=3, color='indianred2') +
  ##   ylab('estimated\ncopy number') +
  ##   xlab('node in collapsed pangenome') +
  ##   facet_grid("coverage\n"~.) +
  ##   scale_y_continuous(breaks=seq(0,10,2)) + 
  ##   theme(legend.position='top', strip.text.y=element_text(angle=0))

  if(fig.scale == 'pangenome'){
    ggp$coverage = ggplot(cov.s, aes(x=node, y=cov.n.2)) + 
      xlab('node in collapsed pangenome')
  } else {
    ggp$coverage = ggplot(cov.s, aes(x=pos, y=cov.n.2)) + 
      xlab('position in chromosome')
  }

  ggp$coverage = ggp$coverage +
    geom_point(alpha=.7) +
    theme_bw() +
    geom_hline(yintercept=4, linetype=1, alpha=.1) +
    geom_hline(yintercept=med.cov.2, linetype=2, alpha=.7) +
    geom_hline(yintercept=fl_cyc_cn, linetype=3, alpha=.7) +
    ylab('estimated\ncopy number') +
    facet_grid("coverage\n"~.) +
    scale_y_continuous(breaks=seq(0,10,2)) + 
    theme(legend.position='top', strip.text.y=element_text(angle=0))

  ## save the x-axis boundaries for later
  if(fig.scale == 'genome') {
    xlims_v = c(xlims_v, cov.s$pos)
  } else {
    xlims_v = c(xlims_v, cov.s$node)
  }

}

##
## haplotypes
##

if(args$viz$val %in% c('diplotype', 'all', 'all_small')){
  ## load stats on pairs of predicted haplotypes
  stats = read.table(args$hstats$val, as.is=TRUE, header=TRUE)
  ## load information for each predicted haplotype
  haps = read.table(args$hpaths$val, as.is=TRUE, header=TRUE)
  ## if same haplotype selected twice, duplicate with a different name
  if(stats$hap1[1] == stats$hap2[1]){
    new_hapn = paste0(stats$hap1[1], '_d')
    stats$hap2[1] = new_hapn
    haps = haps %>% filter(hap==stats$hap1[1]) %>%
      mutate(hap=new_hapn) %>% rbind(haps)
  }
  ## keep info for the best pair only
  haps = haps %>% filter(hap %in% c(stats$hap1[1], stats$hap2[1]))
  
  ## prep predicted haplotype paths
  ggp.haps.df = haps %>% mutate(class=factor(class, c('c1', 'none', 'c2'), c('1', 'both', '2')),
                                ntype=ifelse(class!='both', 'module-specific', 'shared'),
                                hap=factor(hap, levels=unique(hap), labels=c('hap_1', 'hap_2'))) %>% 
    filter(!is.na(class)) %>% arrange(class=='1', class=='2')
  ## add position for each node
  ggp.haps.df = ninfo %>% select(node, pos) %>% merge(ggp.haps.df)
  
  ## how much to vertically shift the points of different classes (for aesthetic purpose)
  path.v.shift.haps = 20
  if(fig.scale == 'pangenome'){
    ggp$haps = ggplot(ggp.haps.df, aes(x=node,
                                       y=ppos + path.v.shift.haps*(as.numeric(class) - 1),
                                       color=class, alpha=ntype)) +
      xlab('node in collapsed pangenome')
  } else {
    ggp$haps = ggplot(ggp.haps.df, aes(x=pos,
                                       y=ppos + path.v.shift.haps*(as.numeric(class) - 1),
                                       color=class, alpha=ntype)) +
      xlab('position in chromosome')
  }

  ## "dot plot" visualization: x=pangenome position, y=haplotype position
  ggp$haps = ggp$haps +
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
           hap=factor(hap, levels=unique(hap), labels=c('hap_1', 'hap_2'))) %>%
    filter(!is.na(ntype))
  path.v.shift.reads = .2

  if(fig.scale == 'pangenome'){
    ggp$haps.s = ggplot(ggp.haps.s.df, aes(x=node,
                                           y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
      xlab('node in collapsed pangenome')
  } else {
    ggp$haps.s = ggplot(ggp.haps.s.df, aes(x=pos,
                                           y=path_part+path.v.shift.reads*(as.numeric(class) - 2))) +
      xlab('position in chromosome')
  }

  ggp$haps.s = ggp$haps.s +
    geom_point(aes(color=class, alpha=ntype)) +
    scale_color_brewer(name='module', palette='Set1') +
    scale_alpha_manual(values=c(.7,.05), name='node') +
    scale_y_continuous(breaks=1:10, expand=c(0,.3)) +
    ylab('path\npart') + 
    theme_bw() + 
    facet_grid(hap~., scales='free', space='free') + 
    theme(legend.position='top',
          axis.text.y=element_text(size=6),
          strip.text.y=element_text(size=4, angle=0)) +
    guides(alpha=guide_legend(ncol=1),
           color=guide_legend(ncol=2),
           shape=guide_legend(ncol=1))

  ## save the x-axis boundaries for later
  if(fig.scale == 'genome') {
    xlims_v = c(xlims_v, ggp.haps.s.df$pos)
  } else {
    xlims_v = c(xlims_v, ggp.haps.s.df$node)
  }

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
xlims = c(min(xlims_v, na.rm=TRUE), max(xlims_v, na.rm=TRUE))

if(fig.scale == 'genome'){
  xlims[2] = xlims[2] + 1000
}

## message('x-axis limits', xlims)
ggp.xlims = xlim(xlims[1], xlims[2])

##
## Genomic element annotation
##

## load gene elements annotation
g.df = read.table(args$genes$val, as.is=TRUE, sep='\t', header=TRUE)

## subset to specific genes if specified in the config
if('genes' %in% names(config)){
  g.df = subset(g.df, gene_name %in% config$genes)
}

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

## adjust to the current (pan)genomic range defined in xlims
if(fig.scale == 'pangenome'){
  g.df = g.df %>% mutate(nstart=ifelse(nstart>xlims[2], xlims[2], nstart),
                         nstart=ifelse(nstart<xlims[1], xlims[1], nstart),
                         nend=ifelse(nend>xlims[2], xlims[2], nend),
                         nend=ifelse(nend<xlims[1], xlims[1], nend)) %>%
    filter(nstart!=nend)
  tss.df = tss.df %>% filter(node<xlims[2], node>xlims[1])
} else {
  g.df = g.df %>% mutate(pstart=ifelse(pstart>xlims[2], xlims[2], pstart),
                         pstart=ifelse(pstart<xlims[1], xlims[1], pstart),
                         pend=ifelse(pend>xlims[2], xlims[2], pend),
                         pend=ifelse(pend<xlims[1], xlims[1], pend)) %>%
    filter(pstart!=pend)
  tss.df = tss.df %>% filter(pos<xlims[2], pos>xlims[1])
}

## start the ggplot object
ggp$genes = g.df %>% filter(type %in% c('exon', 'gene')) %>% 
  ggplot(aes(color=module))

## if variants are also shown, add the vertical lines to help follow their positions
if(!is.null(var.vl)){
  ggp$genes = ggp$genes + var.vl
}

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

if(args$viz$val == 'diplotype'){
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
