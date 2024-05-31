#Taking cluster motifs defined via GIANA on CD154+ CMV or HMGCR-specific TCRs. 
#Projecting those rules onto (a) CD154+ anti-CMV TCRs, (b) CD154+ anti-HMGCR TCRs, (c) HMGCR Myositis Biopsies
#AAG 28 August 2023

setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
library(dplyr)
library(ggplot2)
library(ragg)
library(ggprism)
library(RColorBrewer)
ALPHA=0.05
NAME='ORIG' #For file saving plots. 
source('./_rfunctions/ProcessAdaptiveFile.R')


#Helper Functions----
format_rules <- function(rules){
  rules$seq = gsub('\\?', '\\.', rules$seq) #replace unknown chars w/ * for grep search
  rules$seq = gsub('\\%', '\\.', rules$seq)
  rules$seq = gsub('\\*', '\\.', rules$seq)
  
  rules$trbv = gsub('\\*.*', '', rules$trbv) #omit allele specifiers
  return(rules)
}

patternProjector <- function(data, patterns, includeV = TRUE){
  #split by subject and normalize data. 
  norm = split(data, data$subject)
  norm = lapply(norm, FUN = function(x){
    x$freq = x$count / sum(x$count)
    return(x)
  })
  
  #For each subject, determine total occurence of TCRs meeting any patterns. 
  #Make sure not to double-count a TCR meeting multiple rules. 
  
  #test data
  #patterns = data.frame(seq = c('ASS', 'GE'), trbv = c(NA, NA))
  #x = data.frame(cdr3 = tmp, trbv = 'TRBV3-1', freq = c(0.1,0.2,0.3,0.4), count = c(1,2,3,4))
  
  res = lapply(norm, FUN = function(x){
    meet_criteria = c('')
    for(p in 1:nrow(patterns)){
      cdr3 = grep(patterns$seq[p], x$cdr3)
      if(includeV & !is.na(patterns$trbv[p])){
        trbv = grep(patterns$trbv[p], x$trbv)
        matched = cdr3[which(cdr3 %in% trbv == TRUE)]
      }else{
        matched = cdr3
      }
      meet_criteria = c(meet_criteria, matched)
    }
    meet_criteria = as.numeric(meet_criteria[-1])
    
    if(length(meet_criteria) > 0){
      total_freq = sum(x$freq[unique(meet_criteria)])
      total_count = sum(x$count[unique(meet_criteria)])
    }else{
      total_freq = 0
      total_count = 0
    }
    
    return(data.frame(count = total_count, freq = total_freq))
  })
  
  res = bind_rows(res, .id = 'Subject')
  return(res)
}
GroupComparator <- function(patterns){
  adaptive = lapply(adaptive_groups, patternProjector, patterns = patterns)
  adaptive = bind_rows(adaptive, .id = 'Source')
  adaptive$protocol = 'Adaptive'
  data = adaptive
  return(data)
}
ResultsPlotter <- function(output){
  require(ggplot2)
  require(ggprism)
  require(RColorBrewer)
  
  cols = c(rep('black', 45), brewer.pal(5, 'Dark2'))
  
  gg1 = ggplot(output[output$protocol == 'Adaptive',], aes(x = Source, y = freq)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.3, aes(color = Subject)) + 
    scale_color_manual(values = cols) + 
    ylab('Projection Score') + 
    theme_prism() + 
    theme(legend.position = 'none') 
  
  gg2 = ggplot(output[output$protocol == 'Adaptive' & output$Source != 'Blood' &  output$Source != 'CMV' ,], aes(x = Source, y = freq)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.3) + 
    ylab('Projection Score') + 
    theme_prism()
  
  gg3 = ggplot(output[output$protocol == 'Larman',], aes(x = Source, y = freq)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.3) + 
    ylab('Projection Score') + 
    theme_prism() + 
    coord_flip()
  
  return(list(gg1,gg2,gg3))
}

#Assemble collected sequencing data across all sources-----
load('./objects/combo-hc-myositis.rda')

combined_data  = data.frame(cdr3 = alldata$cdr3_aa, trbv = alldata$v_call, trbj = alldata$j_call, 
                subject = alldata$Sample, source = alldata$Source, count = alldata$duplicate_count, protocol = 'Adaptive')
combined_data$trbv = ProcessAdaptiveVgenes(combined_data$trbv)
periphery = readRDS("./objects/2023-09-17_HMGCR-periphery.rds")

combined_data = combined_data[combined_data$cdr3 != '',]
combined_data = combined_data[-grep('\\*', combined_data$cdr3),] #remove stop codon clonotypes. 
combined_data = combined_data[!is.na(combined_data$cdr3),]

periphery = periphery[periphery$cdr3 != '',]
periphery = periphery[!is.na(periphery$cdr3),]


#Assemble rules to define CMV or HMGCR specific TCRs-----

hmgcr_cluster_results= read.csv("./MotifGroups/2023-09-18_hmgcr-cluster-results-filt_.csv") #GIANA
cmv_cluster_results = read.csv( "./MotifGroups/2023-09-18_cmv-cluster-results-filt.csv") 
hmgcr_cluster_results2 = read.csv("./MotifGroups/2023-11-28_top-giana-clusters-hmgcr-exclusive.csv") #GIANA
load("./objects/2023-09-18_Scrambled-ClusterResults.rda") #Null distributions

#remove common motifs from cmv cluster res prior to projection
hmgcr_cluster_results$motif = paste0(hmgcr_cluster_results$vConsensus, '_', hmgcr_cluster_results$Consensus)
common = hmgcr_cluster_results$motif[which(hmgcr_cluster_results$padj < ALPHA & !(hmgcr_cluster_results$cluster %in% hmgcr_cluster_results2$cluster))]
idx = which((hmgcr_cluster_results$motif %in% common))
hmgcr_scrambles = lapply(hmgcr_scrambles, FUN = function(x){
  x = x[-idx,] #Remove null motif groups which were matched to motifs common with CMV group. 
  return(x)
})

hmgcr_cluster_results = hmgcr_cluster_results2 #Remove motifs common to CMV CD154 group. 
hmgcr_rules = data.frame(seq = hmgcr_cluster_results$Consensus[hmgcr_cluster_results$padj < ALPHA],
                             trbv = hmgcr_cluster_results$vConsensus[hmgcr_cluster_results$padj < ALPHA])     
cmv_cluster_results$motif = paste0(cmv_cluster_results$vConsensus, '_', cmv_cluster_results$Consensus)
cmv_cluster_results = cmv_cluster_results[cmv_cluster_results$padj < ALPHA, ]
cmv_cluster_results_filt = cmv_cluster_results[!(cmv_cluster_results$motif %in% common), ]
cmv_cluster_results_filt = cmv_cluster_results_filt[cmv_cluster_results_filt$padj < ALPHA, ]
idx = which((cmv_cluster_results$motif %in% common))
cmv_scrambles = lapply(cmv_scrambles, FUN = function(x){
  x = x[-idx,] #Remove null motif groups which were matched to motifs common with HMGCR group. 
  return(x)
})


cmv_cluster_results = cmv_cluster_results_filt

idx = match(hmgcr_rules$seq, c$Consensus)
hmgcr_rules_by_clust = split(hmgcr_rules, c$clusters_final[idx])
names(hmgcr_rules_by_clust) = paste0('Group', c(1:8))

hmgcr_scramble_rules = lapply(hmgcr_scrambles, FUN = function(x){
  data.frame(seq = x$Consensus, trbv = x$vConsensus)
})

cmv_rules = data.frame(
            seq = cmv_cluster_results$Consensus[cmv_cluster_results$padj < ALPHA],
            trbv = cmv_cluster_results$vConsensus[cmv_cluster_results$padj < ALPHA])

cmv_scramble_rules = lapply(cmv_scrambles, FUN = function(x){
  data.frame(seq = x$Consensus, trbv = x$vConsensus)
})
#Code----

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
# #Before projecting, add Biopsy and Periphery from autoFEST subject.
# fnames = paste0('HB_Baseline_', c(1:4), '.txt')
# full_paths = paste0('./data/ET_myositis/2023_10_ET-Autofest-data/', fnames)
# m = lapply(full_paths, read.delim)
# m = bind_rows(m) #I don't care which file it came from. 
# 
# m6_periphery  = data.frame(
#   cdr3 = m$cdr3aa, 
#   trbv = m$v, 
#   trbj = m$j, 
#   subject = 'HB', 
#   source = 'Periphery',
#   count = m$count,
#   protocol = 'KS-core'
# )
# 
# m = read.delim('./data/ET_myositis/2023_10_ET-Autofest-data/R23-923.txt')
# m6_biopsy = data.frame(
#   cdr3 = m$cdr3aa, 
#   trbv = m$v, 
#   trbj = m$j, 
#   subject = 'HB', 
#   source = 'Biopsy',
#   count = m$count,
#   protocol = 'KS-core'
# )
# 
# adaptive_groups$Biopsy = rbind(adaptive_groups$Biopsy, m6_biopsy)
# adaptive_groups$Periphery = rbind(adaptive_groups$Periphery, m6_periphery)

#Now proceed with HMGCR Projection------
rule_sets = hmgcr_rules_by_clust

rule_sets = lapply(rule_sets, format_rules)#Format rule sets for grep search. 
results = lapply(rule_sets, GroupComparator) 

results = bind_rows(results, .id = 'Motif Group')
results$Group = substr(results$`Motif Group`, 6, nchar(results$`Motif Group`))
#Make Stacked Bar Plot of Projections by Group. 
hmgcr = results[results$Source == 'Blood',]
biopsy = results[results$Source == 'Biopsy',]
cmv = results[results$Source == 'CMV',]
hc = results[results$Source == 'HC_PBMC',]
#cols = brewer.pal(9, 'Set3')[-2]
#cols = brewer.pal(8, 'Dark2')[c(1:6, 8, 7)]
cols = c("#66C2A5" ,"#8DA0CB" ,"#E78AC3", "#A6D854", "#FC8D62" , "#FFD92F",  "#B3B3B3", "#E5C494") 
#cols = brewer.pal(8, 'Set2')[c(1:3,5:6, 8, 7,4)] 
gg1 = ggplot(hmgcr, aes(x = Subject, y = freq * 100, fill = Group)) + 
         geom_col() + 
          theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'CD154+ anti-HMGCR TCRs') + 
  theme(axis.text.x = element_text(size = 8))

gg2 = ggplot(biopsy, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'Muscle Biopsy TCRs') + 
  theme(axis.text.x = element_text(size = 8))

gg3 = ggplot(cmv, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'CD154+ anti-CMV TCRs') + 
  theme(axis.text.x = element_text(size = 8))

gg4 = ggplot(hc, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'HC Peripheral TCRs') + 
  theme(axis.text.x = element_text(8)) + 
  ylim(0,1.5)
  
plots = list(hmgcr_hclust = gg1, biopsy_hclust = gg2, cmv_hclust= gg3, 
             hc_hclust = gg4)
plots_hmgcr = plots
for(p in c(1:length(plots))){

 
  fnames = paste0('./', Sys.Date(),  '_', names(plots)[p], c(1:3))
  agg_png(filename = fnames[1], width = 1600, height = 1200, units = 'px', scaling = 3)
  plot(plots[[p]])
  dev.off()
  
  agg_png(filename = fnames[2], width = 800, height = 1200, units = 'px', scaling = 3)
  plot(plots[[p]])
  dev.off()
  
  agg_png(filename = fnames[3], width = 1200, height = 800, units = 'px', scaling = 3)
  plot(plots[[p]])
  dev.off()
}

save(results, plots_hmgcr, file = paste0('./', Sys.Date(), '_ProjectionResults-exclusive.rda'))

#Now proceed with CMV Projection------
rule_sets = c(list(cmv = cmv_rules), cmv_scramble_rules)

rule_sets = lapply(rule_sets, format_rules)#Format rule sets for grep search. 
results = lapply(rule_sets, GroupComparator) 
results_orig = results

#produce mean of 10 scambles 
cmv_scramble_res = results[2:11]

tmp = bind_rows(cmv_scramble_res, .id = 'Trial')
spl = split(tmp, paste0(tmp$Source, '_', tmp$Subject))
spl = lapply(spl, FUN = function(x){
  data.frame(Source = x$Source[1], Subject = x$Subject[1], 
             freq = mean(x$freq), count  = mean(x$count), protocol = x$protocol[1])
})
tmp = bind_rows(spl)
tmp$Group = 'Scramble'
#tapply(tmp$freq, paste0(tmp$Subject, FUN = mean)



results = results$cmv
results$Group = 'CMV'
results = bind_rows(results, tmp )

tmp = results[results$Source %in% c('Biopsy', 'HC_PBMC'),]
tmp$x = paste0(tmp$Source,'_', tmp$Group)
tmp$x = factor(tmp$x, levels = c('HC_PBMC_Scramble', 'HC_PBMC_CMV', 
               'Biopsy_Scramble', 'Biopsy_CMV'))
ggplot(tmp, aes(x = x, y = freq * 100, color = Source)) + 
  geom_boxplot() + 
  theme_prism() + 
  theme(axis.text.x = element_text(angle = 90))
#Make Stacked Bar Plot of Projections by Group. 
hmgcr = results[results$Source == 'Blood',]
biopsy = results[results$Source == 'Biopsy',]
cmv = results[results$Source == 'CMV',]
hc = results[results$Source == 'HC_PBMC',]
per = results[results$Source == 'Periphery',]
cols = brewer.pal(9, 'Set3')[-2]

gg1 = ggplot(hmgcr, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'CD154+ anti-HMGCR TCRs') + 
  theme(axis.text.x = element_text(size = 11))

gg2 = ggplot(biopsy, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'Muscle Biopsy TCRs') + 
  theme(axis.text.x = element_text(size = 11))

gg3 = ggplot(cmv, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'CD154+ anti-CMV TCRs') + 
  theme(axis.text.x = element_text(size = 11))

gg4 = ggplot(hc, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'HC Peripheral TCRs') + 
  theme(axis.text.x = element_blank()) + 
  ylim(0,1.5)

gg5 = ggplot(per, aes(x = Subject, y = freq * 100, fill = Group)) + 
  geom_col() + 
  theme_prism() + 
  scale_fill_manual(values = cols) + 
  ylab('Projection Score') + 
  labs(title = 'Myositis Peripheral TCRs') 

plots = list(hmgcr_hclust = gg1, biopsy_hclust = gg2, cmv_hclust= gg3, 
             hc_hclust = gg4, m_hclust = gg5)
plots_cmv = plots
for(p in c(1:length(plots))){
  
  
  fnames = paste0('./', Sys.Date(),  '_', names(plots)[p], c(1:3))
  agg_png(filename = fnames[1], width = 1600, height = 1200, units = 'px', scaling = 3)
  plot(plots[[p]][[1]])
  dev.off()
  
  agg_png(filename = fnames[2], width = 800, height = 1200, units = 'px', scaling = 3)
  plot(plots[[p]][[2]])
  dev.off()
  
  agg_png(filename = fnames[3], width = 1200, height = 800, units = 'px', scaling = 3)
  plot(plots[[p]][[3]])
  dev.off()
}

save(results_orig,  results, plots_cmv, file = paste0('./', Sys.Date(), '_ProjectionResults-CMV-exclusive.rda'))

