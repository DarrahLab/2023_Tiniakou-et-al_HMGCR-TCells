setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
load('./objects/combo-hc-myositis.rda')
source('./_rfunctions/ProcessAdaptiveFile.R')
library(rdi)
library(umap)
set.seed(99)

shannon_equitability <- function(freqs){
  freqs = freqs/sum(freqs)
  tmp = lapply(freqs, FUN = function(x){
    return(-x * log(x))
  })
  shannon = 1 - sum(unlist(tmp)) / log(length(freqs))
  return(shannon)
}
shannon_subsampler <- function(freqs, seqs = NA, N = 1e3, I = 1){
  #Input sequences at their native frequency: more common sequences 
  #should appear more often. Otehrwise supply freqs. N is the number of seqs to subsample 
  # I is the number of iterations. 
  i = 0
  shannon_eqs = vector(mode = 'double', length = I)
  if(all(is.na(seqs))) seqs = paste0('Seq', c(1:length(freqs)))
  
  while(i < I){
    i = i + 1
    samp = sample(seqs, N, prob = freqs, replace = TRUE)
    samp_freqs = as.numeric(table(samp) / sum(table(samp)))
    shannon_eqs[i] = shannon_equitability(samp_freqs)
  }
  
  return(mean(shannon_eqs))
}
#Assemble collected sequencing data across all sources-----

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
 data$productive = data$cdr3 != '' & !is.na(data$cdr3) & !(c(1:nrow(data)) %in% grep('\\*', data$cdr3))
tab1 = table(data$productive, data$Group)
tab2 = apply(tab1, 2, FUN = function(x) round(x/sum(x) * 100,1))
#Only consider gene 
data$trbv = gsub("\\*.*","",data$trbv)
data$uid = paste0(data$trbv, '_', data$cdr3)

spl = split(data, data$Group)
x = spl$Periphery
y = split(x, x$subject)
y = lapply(y, FUN = function(q){
  q = q[sample(c(1:nrow(q)), 1e4),]
  return(q)})
x = bind_rows(y)
spl$Periphery_subset = x
spl= lapply(spl, FUN = function(x){
  
  y = split(x, x$subject)
  y = lapply(y, FUN = function(q){
    counts = tapply(q, q$uid, FUN = function(z){sum(z$count)})
    se = shannon_subsampler(counts)
    return(se)
  })
 return(unlist(y))
})
res = vector(mode = 'list', length = length(spl))
for(s in 1:length(spl)){
  res[[s]] = data.frame(Source = names(spl)[s], equitability = as.numeric(spl[[s]]))
}

res = bind_rows(res)

plot1 = ggplot(res, aes(x = Source, y = equitability))  + 
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.3) + 
  #geom_violin() + 
  geom_jitter(position = position_jitter( width = 0.25), size = 0.75, alpha = 1, aes(color = Source))  

#res_orig = res
res = vector(mode = 'list', length = length(spl))
for(s in 1:length(spl)){
  res[[s]] = data.frame(Source = names(spl)[s], equitability = as.numeric(spl[[s]]))
}

res = bind_rows(res)

res  = res[res$Source %in% c('Biopsy', 'HC_PBMC', 'Periphery_subset'),]
res$Source[res$Source == 'HC_PBMC'] = 'Blood HC'
res$Source[res$Source == 'Periphery_subset'] = 'Blood IMNM'
res$Source = factor(res$Source, levels = c('Blood HC', 'Blood IMNM', 'Biopsy'))

cols = brewer.pal(3, 'Set2')[c(1,3,2)]
plot2 = ggplot(res, aes(x = Source, y = equitability))  + 
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.3) + 
  #geom_violin() + 
  geom_jitter(position = position_jitter( width = 0.25), size = 0.75, alpha = 1, 
              aes(color = Source))  + 
  scale_color_manual(values = cols) + 
  theme_prism() + 
  guides(colour = guide_legend(override.aes = list(size=3)))
  
plot2 

#This is the final plot to be included in supplement. 


require(ggpubr)
plot3 = ggplot(res, aes(x = Source, y = equitability))  + 
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.3) + 
  #geom_violin() + 
  geom_jitter(position = position_jitter( width = 0.25), size = 1, alpha = 1, 
              aes(color = Source))  + 
  scale_color_manual(values = c('#12405A', '#E9AA74', '#E8825D'))+  
  theme_prism() + 
  theme(axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  ylab('Shannon Equitability') + 
  stat_compare_means(method = 't.test', 
                     comparisons = list(
                       
                       c('Biopsy', 'Blood IMNM'), 
                       c('Blood HC', 'Blood IMNM'), 
                       c('Biopsy', 'Blood HC'))) 
plot3

write.csv(res, file = paste0('./code-outputs/',Sys.Date(), '_shannon-eq.csv'))
