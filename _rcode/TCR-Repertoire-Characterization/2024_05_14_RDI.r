#Obtain RDI for Eleni myositis samples as well as Lyme cohort healthy controls. 
#AAG 16 February 2023
setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells')
load('./objects/combo-hc-myositis.rda')
source('./_rfunctions/ProcessAdaptiveFile.R')
library(rdi)
library(umap)
set.seed(99)
#Assemble collected sequencing data across all sources-----
load('./objects/combo-hc-myositis.rda')


combined_data = data.frame(cdr3 = alldata$cdr3_aa, trbv = alldata$v_call, trbj = alldata$j_call, 
                subject = alldata$Sample, source = alldata$Source, count = alldata$duplicate_count, protocol = 'Adaptive')


combined_data$trbv = ProcessAdaptiveVgenes(combined_data$trbv)
periphery = readRDS("./objects/2023-09-17_HMGCR-periphery.rds")

spl = split(combined_data, combined_data$source)
adaptive_groups = spl[which(names(spl) %in% unique(combined_data$source[combined_data$protocol == 'Adaptive']))]
adaptive_groups$Periphery = data.frame(
  cdr3 = periphery$cdr3, 
  trbv = periphery$trbv, 
  trbj = NA, 
  subject = periphery$sample, 
  source = 'Periphery', 
  count = 1, 
  protocol = 'Adaptive'
)

data = bind_rows(adaptive_groups, .id = 'Group')
data = data[data$cdr3 != '',]
data = data[-grep('\\*', data$cdr3),] #remove stop codon clonotypes. 
data = data[!is.na(data$cdr3),]
#Only consider gene 
data$trbv = gsub("\\*.*","",data$trbv)
data$ID = paste0(data$subject, '_', data$Group)


genes = data.frame(rep(data$trbv, 
                       data$count))
seqAnnot = rep(data$ID, data$count)
cts = calcVDJcounts(genes,seqAnnot) 

cts = cts[rownames(cts) != 'NA-1',]
rdi <- calcRDI(cts, units = "pct")

#Computing pairwise distances, are Biopsies closer to eachother than by chance? #####
#I mean obviously yes but prove it. 
d = as.matrix(rdi)
bidx = grep('Biopsy', rownames(d))
hcidx = grep('HC_PBMC', rownames(d))
hmgcridx = grep('Blood', rownames(d))
cmvidx = grep('CMV', rownames(d))

groups = c('Biopsy', 'Periphery', 'HC')
#groups = c('Biopsy', 'HC')
sim_dfs = lapply(groups, FUN = function(x){
  idx = grep(x, rownames(d))
  d_sub = d[idx, idx]
  d_sub[upper.tri(d_sub, diag = TRUE)] = NA
  df = reshape2::melt(d_sub)
  df = df[-which(is.na(df$value)),]
  return(df)
})
names(sim_dfs) = groups
sims = bind_rows(sim_dfs, .id = 'Group')

save(sims, file = paste0('./code-outputs/',Sys.Date(), '_rdi.rda'))
