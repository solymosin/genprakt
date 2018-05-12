# R

library(VariantAnnotation)
library(GenomicFeatures)

options(stringsAsFactors=F)

VPhaser2.parser = function(fajl, seqnames){
  raw.txt = readLines(fajl)
  raw.txt = raw.txt[-grep('# ---', raw.txt)]
  raw.txt = raw.txt[-grep('# Summary:', raw.txt)]
  fejlec = raw.txt[grep('# Ref_Pos', raw.txt)]
  fejlec = strsplit(gsub('# ', '', fejlec), '\t')[[1]]
  raw.txt = raw.txt[-grep('# Ref_Pos', raw.txt)]

  i = 1
  sor = strsplit(raw.txt[i], '\t')[[1]]
  sorok = c(sor[1:6],paste(sor[7:length(sor)], collapse=' '))
  for(i in 2:length(raw.txt)){
      sor = strsplit(raw.txt[i], '\t')[[1]]
      sorok = rbind(sorok, c(sor[1:6],paste(sor[7:length(sor)], collapse=' ')))
  }

  rownames(sorok)=1:dim(sorok)[1]
  tab = data.frame(sorok)
  colnames(tab) = fejlec
  tab$Ref_Pos = as.numeric(tab$Ref_Pos)
  tab$Strd_bias_pval = as.numeric(tab$Strd_bias_pval)
  tab$Var_perc = as.numeric(tab$Var_perc)
  
  tab$seqnames = seqnames
  tab$start = tab$Ref_Pos
  tab$end = tab$Ref_Pos
  tab$strand = '*'
  tab$REF = tab$Cons
  tab$ALT = tab$Var
  oszlopok = c('seqnames', 'start', 'end', 'strand', 'REF', 'ALT', 'Strd_bias_pval', 'Var_perc', 'SNP_or_LP_Profile', 'Type')
  gr = GRanges(tab[,oszlopok])
  return(gr)  
}



