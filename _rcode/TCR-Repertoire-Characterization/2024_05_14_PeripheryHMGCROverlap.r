#Goal: Identify the percentage of peripheral blood and biopsy repertories which 
# are comprised of CD4+ CD14+ anti-HMGCR TCRs. 
#AAG 05/2024

setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
load('./objects/combo-hc-myositis.rda')
source('./_rfunctions/ProcessAdaptiveFile.R')
library(rdi)
library(umap)
library(dplyr)
library(parallel) #clone tables computation. 
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggpubr)
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
hmgcr = hmgcr[!(hmgcr$uid %in% in_common),]

#data = data[!(data$uid %in% in_common),] #remove non ag specific clonotypes
data = data[data$cdr3 != '',]
data = data[-grep('\\*', data$cdr3),] #remove stop codon clonotypes. 
data = data[!is.na(data$cdr3),]
hmgcr = data[data$Group == 'Blood',]
hmgcr = hmgcr[!(hmgcr$uid %in% in_common),]
biopsy = data[data$source == 'Biopsy',]
spl1 = split(biopsy, biopsy$subject)
per = data[data$source == 'Periphery',]
spl2 = split(per, per$subject)

spl1 = lapply(spl1, FUN = function(x){
  x = x[match(unique(x$uid), x$uid),]
  sum(x$uid %in% hmgcr$uid) / nrow(x) * 100
})

# spl1 = lapply(spl1, FUN = function(x){
#   sum(x$uid %in% hmgcr$uid) #/ nrow(x) * 100
# })


spl2 = lapply(spl2, FUN = function(x){
  x = x[match(unique(x$uid), x$uid),]
  sum(x$uid %in% hmgcr$uid) / nrow(x) * 100
})

spl2 = t(bind_rows(spl2, .id = 'Subject'))
spl1 = t( bind_rows(spl1, .id = 'Subject'))

df = data.frame(
  Subject = as.character(c(12211,13038,17137,18088,18124,
                           12211,13038,17137,18124)), 
  Source = c(rep('Biopsy', 5), rep('Periphery', 4)), 
  Frequency = c(0.7936508, 0.2818035, 0, 0.2556455, 0.2619515, 
                0.033662306,0.054185857, 0.006884885,  0.003213058)
)
df$Source = factor(df$Source, levels = c('Periphery', 'Biopsy'))
ggplot(df, aes(x = Source, y = Frequency)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(color = Subject)) +
  scale_color_brewer(palette = 'Dark2') + 
  stat_compare_means(method = 't.test') + 
  theme_prism()
