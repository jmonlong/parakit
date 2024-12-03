library(dplyr)

##
## ClinVar

## read ClinVar for specific columns
df = read.table('variant_summary.txt.gz', sep='\t', as.is=TRUE, header=TRUE, comment.char='', quote='', nrows=1)
colc = rep('NULL', ncol(df))
cols.to.keep = c('Chromosome', 'Start', 'Stop', 'X.AlleleID',
                 'ReferenceAlleleVCF', 'AlternateAlleleVCF',
                 'GeneSymbol',
                 'Name', 'ClinicalSignificance', 'Assembly')
colc[which(colnames(df) %in% cols.to.keep)] = 'character'

df = read.table('variant_summary.txt.gz', sep='\t', colClasses=colc, header=TRUE, comment.char='', quote='')

## keep variants with either the word "Pathogenic" or "Conflicting" in the
## clinical significance column, and the information based on the GRCh38
## reference genome.
df.p = df %>% filter(sapply(strsplit(GeneSymbol, ';'), function(gg) any(gg == 'CYP21A2')), 
                     grepl("Conflicting", ClinicalSignificance) |
                     grepl("Pathogenic", ClinicalSignificance),
                     Assembly=='GRCh38') %>%
  mutate(nuc.change=gsub('.*:c.([^ ]*).*', '\\1', Name),
         prot.change=gsub('.*:.* \\(p.(.*)\\)', '\\1', Name),
         prot.change=ifelse(grepl('p.', Name), prot.change, NA)) %>%
  rename(chr=Chromosome, start=Start, end=Stop, id=X.AlleleID,
         ref=ReferenceAlleleVCF, alt=AlternateAlleleVCF,
         clnsig=ClinicalSignificance) %>% 
  select(chr, start, end, ref, alt, id, nuc.change, prot.change, clnsig)
sample_n(df.p, 10)

## read the date when the clinvar information was downloaded
cv.date = scan('variant_summary.txt.date', '', quiet=TRUE)

write.table(df.p, row.names=FALSE, quote=FALSE, sep='\t',
            file=paste0('CYP21A2.pathogenic.variant_summary.', cv.date, '.tsv'))

##
## gene annotation

## download GENCODE if necessary
if (!file.exists('gencode.v43.basic.annotation.gtf.gz')){
  download.file('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz', 'gencode.v43.basic.annotation.gtf.gz')
}

## read GENCODE
df = read.table('gencode.v43.basic.annotation.gtf.gz', as.is=TRUE, sep='\t')
df = df[,c(1,4,5,3,7,9,2)]
colnames(df) = c('chr', 'start', 'end', 'type', 'strand', 'info', 'source')

## parse gene name
df = df %>% mutate(gene_name=gsub('.*gene_name ([^;]+).*', '\\1', info))

## genes to keep
## either around a position
## df %>% filter(chr=='chr6', abs(start-32040869) < 6e4, type=='gene')
## or specify gene names directly
sel.genes = c('C4A', 'CYP21A1P', 'TNXA', 'C4B', 'CYP21A2', 'TNXB')

## keep information about the UTR/exons and coding sequence
sel.types = c('gene', 'exon', 'utr', 'CDS')
df.g = df %>% filter(gene_name %in% sel.genes, type %in% sel.types) %>%
  mutate(gene_type=gsub('.*gene_type ([^;]+).*', '\\1', info))

## write annotation file
df.g %>% select(-info, -source) %>% unique %>%
  write.table('CYP21A2.gencodev43.nearby_genes.tsv', sep='\t', row.names=FALSE, quote=FALSE)
