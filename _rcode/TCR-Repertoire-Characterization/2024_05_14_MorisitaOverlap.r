#Obtain Morisita overlap and the size of intersection between each set of repertories: 
# 45 HC repertoires, HMGCR periphery, HMGCR muscle biopsies, CD4+CD154+ populations. 
#AAG 16 February 2023
setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
source('./_rfunctions/ProcessAdaptiveFile.R')
library(rdi)
library(umap)
library(dplyr)
library(parallel) #clone tables computation. 
library(ggplot2)
library(ggprism)
library(RColorBrewer)
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

source('./_rfunctions/ProcessAdaptiveFile.R')
data$trbv = ProcessAdaptiveVgenes(data$trbv)

#Only consider gene 
data$trbv = gsub("\\*.*","",data$trbv)
data$ID = paste0(data$subject, '_', data$Group)

#Remove TCRs in common to CD154+CD4+ anti-CMV and anti-HMGCR populations and nonproductives####
data$uid = paste0(data$trbv, '_', data$cdr3)
hmgcr = data[data$Group == 'Blood',]
cmv = data[data$Group == 'CMV',]
in_common = hmgcr$uid[hmgcr$uid %in% cmv$uid]
common_clonotypes = in_common[setdiff(c(1:95), grep('NA', in_common))]
all = data
tab_all = table(all$subject[all$Group %in% c('Blood', 'Biopsy', 'Periphery', 'CMV')], 
      all$Group[all$Group %in% c('Blood', 'Biopsy', 'Periphery', 'CMV')])
require(xlsx)
write.xlsx(as.data.frame.matrix(tab_all), sheetName = 'AllTCRs', 
           file = './code-outputs/Myositis-TCR-Yields.xlsx')
data = data[data$cdr3 != '',]
data = data[-grep('\\*', data$cdr3),] #remove stop codon clonotypes. 
data = data[!is.na(data$cdr3),]
prod = data
tab_prod = table(prod$subject[prod$Group %in% c('Blood', 'Biopsy', 'Periphery', 'CMV')], 
      prod$Group[prod$Group %in% c('Blood', 'Biopsy', 'Periphery', 'CMV')])
write.xlsx(as.data.frame.matrix(tab_prod), sheetName = 'ProductiveTCRs', 
           file = './code-outputs/Myositis-TCR-Yields.xlsx',append = TRUE)
data = data[!(data$uid %in% in_common),] #remove non ag specific clonotypes
###Run Similarity Stuff######
data$Sample = paste0(data$subject, '_', data$source)
yields = tapply(data, data$Sample, FUN = function(x){return(sum(x$count))})
data$clone_id = paste0(data$trbv, '_', data$cdr3)


#number of subjects per clonotype#####
subset = prod[prod$Group == 'Biopsy',]
tab = table(subset$uid, subset$subject) > 0
nSubjects = rowSums(tab)
tab2 = table(nSubjects)
sdf = data.frame(nSubjects = c(1:5), nClonotypes = as.numeric(tab2))
write.csv(sdf, row.names = FALSE, quote = FALSE, file = './code-outputs/Biopsy-public-count.csv')
ggplot(sdf, aes(x = nSubjects, y = nClonotypes)) + 
  geom_col() + 
  theme_prism() + 
  scale_y_log10()
####pre-computing specific things####
#Make clone tables of common clonotypes between all subjects
spl = split(data, data$Sample)
clone_tables = mclapply(spl, mc.cores = 8, FUN = function(x){
  tab = tapply(x, x$clone_id, FUN = function(y){sum(y$count)})
  return(tab)
})
simpson <- function(x){
  #x in vector or table of unique clonotype counts. 
  p = x/sum(x)
  s = sum(p^2)
  return(s)
}
simpson_index = mclapply(clone_tables, mc.cores = 8, FUN = simpson)
#Run overlaps######
#If I were smart I would have parallelized this, but oh well. 
# Overlap Coefficient: a normalised measure of overlap similarity. 
#It is defined as the size of the intersection divided by the smaller of 
#the size of the two sets.
overlap_coeff = matrix(data = NA, nrow = length(clone_tables), 
                        ncol = length(clone_tables))

#NSeqs: The number of clonotypes in common between two sets. 
nseq_matrix = matrix(data = NA, nrow = length(clone_tables), 
                        ncol = length(clone_tables))
#CSeqs: The number of total cells with a sequence common between two sets. 
cseq_matrix = matrix(data = NA, nrow = length(clone_tables), 
                     ncol = length(clone_tables))
#Morisita Overlap: 
morisita_matrix = matrix(data = NA, nrow = length(clone_tables), 
                     ncol = length(clone_tables))

for(i in 1:length(clone_tables)){
  for(j in 1:length(clone_tables)){
    all_clones = unique(c(names(clone_tables[[i]]), names(clone_tables[[j]])))
    data = data.frame(
      clone = c(names(clone_tables[[i]]), names(clone_tables[[j]])), 
      count = c(as.numeric(clone_tables[[i]]), as.numeric(clone_tables[[j]])), 
      sample = c(rep(names(clone_tables)[i], length(clone_tables[[i]])), 
                     rep(names(clone_tables)[j], length(clone_tables[[j]])  )) 
    )
  
    tab = table(data$clone, data$sample)
    logi = tab > 0
    
    idx = which(rowSums(logi ) >1)
    
    num = 2 * sum(unlist(lapply(idx, FUN = function(x){tab[x,1] * tab[x,2]})))
    den = (simpson_index[[i]] + simpson_index[[j]]) * sum(clone_tables[[i]]) * sum(clone_tables[[j]])
    
    overlap_coeff[i,j] = sum(tab[idx,]) / sum(clone_tables[[i]])
    nseq_matrix[i,j] = length(idx)
    cseq_matrix[i,j] = sum(tab[idx,1])
    morisita_matrix[i,j] = num/den
  }
}

rownames(morisita_matrix) = names(clone_tables); colnames(morisita_matrix) = names(clone_tables)
rownames(nseq_matrix) = names(clone_tables); colnames(nseq_matrix) = names(clone_tables)
rownames(cseq_matrix) = names(clone_tables); colnames(cseq_matrix) = names(clone_tables)
rownames(overlap_coeff) = names(clone_tables); colnames(overlap_coeff) = names(clone_tables)

overlap_coeff[upper.tri(overlap_coeff)] = NA
nseq_matrix[upper.tri(nseq_matrix)] = NA
morisita_matrix[upper.tri(morisita_matrix)] = NA

ovdf = reshape2::melt(overlap_coeff)
nsdf = reshape2::melt(nseq_matrix)
msdf = reshape2::melt(morisita_matrix)

ovdf$metric = 'Overlap Coefficient'
nsdf$metric = 'Overlapping Clonotypes'
msdf$metric = 'Morisita Overlap'

metric = rbind(ovdf, nsdf, msdf)

torm = which(metric$Var1 == metric$Var2)
metric = metric[-torm,]
# torm = c(grep('NA', metric$Var1), grep('NA', metric$Var2))
# metric=metric[-torm,]
metric$Var1 = as.character(metric$Var1); metric$Var2 = as.character(metric$Var2)
metric = metric[which(!is.na(metric$value)),]

metric$Group1 = NA; metric$Group2 = NA; 

metric$Group1[grep('HC_PBMC', metric$Var1)] = 'HC'
metric$Group1[grep('Periphery', metric$Var1)] = 'Periphery'
metric$Group1[grep('Biopsy', metric$Var1)] = 'Biopsy'
metric$Group1[grep('Blood', metric$Var1)] = 'HMGCR'
metric$Group1[grep('CMV', metric$Var1)] = 'CMV'

metric$Group2[grep('HC_PBMC', metric$Var2)] = 'HC'
metric$Group2[grep('Periphery', metric$Var2)] = 'Periphery'
metric$Group2[grep('Biopsy', metric$Var2)] = 'Biopsy'
metric$Group2[grep('Blood', metric$Var2)] = 'HMGCR'
metric$Group2[grep('CMV', metric$Var2)] = 'CMV'

metric$Comparison = paste0(metric$Group1, '_', metric$Group2)

metric$Comparison[metric$Comparison == 'CMV_Biopsy'] = 'Biopsy_CMV'
metric$Comparison[metric$Comparison == 'Periphery_Biopsy'] = 'Biopsy_Periphery'
metric$Comparison[metric$Comparison == 'HMGCR_Biopsy'] = 'Biopsy_HMGCR'
metric$Comparison[metric$Comparison == 'CMV_Periphery'] = 'Periphery_CMV'
metric$Comparison[metric$Comparison == 'HMGCR_Periphery'] = 'Periphery_HMGCR'
metric$Comparison[metric$Comparison == 'HMGCR_CMV'] = 'CMV_HMGCR'

save(metric, morisita_matrix, nseq_matrix, cseq_matrix, overlap_coeff, 
     file = paste0('./code-outputs/', Sys.Date(), '_similarity-metrics.rda'))
#Two public CD154+ clonotypes found in those of another subject. 
tmp = metric[metric$Comparison == 'HMGCR_HMGCR' & metric$metric == 'Overlapping Clonotypes',]
tmp = metric[metric$Comparison == 'Biopsy_HMGCR' & metric$metric == 'Overlapping Clonotypes',]
tmp = tmp[,1:3]
tmp = metric[metric$Comparison == 'Biopsy_Biopsy' & metric$metric == 'Overlapping Clonotypes',]











