#Taking cluster motifs defined via GIANA on CD154+ CMV or HMGCR-specific TCRs. 
#Projecting those rules onto (a) CD154+ anti-CMV TCRs, (b) CD154+ anti-HMGCR TCRs, (c) HMGCR Myositis Biopsies, 
#(d) Mammen Data. 
#Difference from before is that I am excluding 6 motifs in common w/CMV. 
#AAG 29 November 2023 

setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
library(dplyr)
library(ggplot2)
library(ragg)
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
consolidateGLIPHResults <- function(gliph){
  require(dplyr)
  spl = split(gliph, gliph$index)
  spl = lapply(spl, FUN = function(x){
    if(length(grep('global', x$type[1]))){
      df = data.frame(seq = gsub('.*\\global-', '', x$type[1]), 
                      trbv = NA, 
                      pvalue = x$final_score[1], 
                      type = 'global'
      )
      
    }else if(length(grep('motif', x$type[1]))){
      local_motif = gsub('motif-', '', x$type[1])
      #Determine positions of local motif in cluster set. 
      positions = unlist(lapply(x$TcRb, FUN = function(y){
        tmp = gregexpr(pattern = local_motif, y, ignore.case = TRUE)
        return(tmp[[1]][1])
      }))
      max_pos = max(positions)
      min_pos = min(positions)
      
      patterns = unlist(lapply(c(min_pos:max_pos) - 1, FUN = function(n){
        paste0( paste0(rep(c('.'), n), collapse = ''), local_motif)
      }))
      
      df = data.frame(seq = patterns, 
                      trbv = NA, 
                      pvalue = x$final_score[1], 
                      type = 'local'
      )
    }
    
    vtab = table(x$V)
    ratio = max(vtab)/sum(vtab)[1]
    if(ratio > 0.5){
      df$trbv = names(vtab)[which(vtab==max(vtab))]
    }
    return(df)
  })
  
  out = bind_rows(spl)
  return(out)
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
  frseq = lapply(larman_groups, patternProjector, patterns = patterns)
  
  adaptive = bind_rows(adaptive, .id = 'Source')
  frseq = bind_rows(frseq, .id = 'Source')
  adaptive$protocol = 'Adaptive'
  frseq$protocol = 'Larman'
  
  data = rbind(adaptive, frseq)
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


combined_data = data.frame(cdr3 = alldata$cdr3_aa, trbv = alldata$v_call, trbj = alldata$j_call, 
                subject = alldata$Sample, source = alldata$Source, count = alldata$duplicate_count, protocol = 'Adaptive')

combined_data$trbv = ProcessAdaptiveVgenes(combined_data$trbv)
periphery = readRDS("./objects/ET-tcr/2023-09-17_HMGCR-periphery.rds")

combined_data = combined_data[combined_data$cdr3 != '',]
combined_data = combined_data[-grep('\\*', combined_data$cdr3),] #remove stop codon clonotypes. 
combined_data = combined_data[!is.na(combined_data$cdr3),]

periphery = periphery[periphery$cdr3 != '',]
periphery = periphery[!is.na(periphery$cdr3),]
#Assemble rules to define CMV or HMGCR specific TCRs-----

#hmgcr_cluster_results = read.csv("./results/ET_myositis/2023-09-18_hmgcr-cluster-results_ORIG.csv") #GIANA
#cmv_cluster_results = read.csv( "./results/ET_myositis/2023-09-18_cmv-cluster-results_ORIG.csv") #GIANA

#hmgcr_cluster_results = read.csv("./results/ET_myositis/2023-09-18_hmgcr-cluster-results-filt_.csv") #GIANA
#cmv_cluster_results = read.csv( "./results/ET_myositis/2023-09-18_cmv-cluster-results-filt.csv") #GIANA
load("./objects/2023-09-18_Scrambled-ClusterResults.rda") #Null distributions
hmgcr_cluster_results = read.csv("./MotifGroups/2023-11-28_top-giana-clusters-hmgcr-exclusive.csv")
 
hmgcr_rules = data.frame(seq = hmgcr_cluster_results$Consensus[hmgcr_cluster_results$padj < ALPHA],
                             trbv = hmgcr_cluster_results$vConsensus[hmgcr_cluster_results$padj < ALPHA])         

hmgcr_scramble_rules = lapply(hmgcr_scrambles, FUN = function(x){
  data.frame(seq = x$Consensus, trbv = x$vConsensus)
})

# cmv_scramble_rules = lapply(cmv_scrambles, FUN = function(x){
#   data.frame(seq = x$Consensus, trbv = x$vConsensus)
# })
common = c(4,6,10,22,23,24) #These motifs in scrambles should be excluded. 

hmgcr_scramble_rules = lapply(hmgcr_scramble_rules, FUN = function(x){
  x = x[-common,]
  return(x)
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
larman_groups = spl[which(names(spl) %in% unique(combined_data$source[combined_data$protocol == 'Larman']))]

rule_sets = list(hmgcr_giana = hmgcr_rules)# ,cmv_giana = cmv_rules)


rule_sets = lapply(rule_sets, format_rules)#Format rule sets for grep search. 
results = lapply(rule_sets, GroupComparator) 
plots = lapply(results, ResultsPlotter)


##Do this next. 
concatScramble <- function(scramble_res){
  Ct = matrix(data = NA, nrow = nrow(scramble_res[[1]]), ncol = length(scramble_res))
  Fr = matrix(data = NA, nrow = nrow(scramble_res[[1]]), ncol = length(scramble_res))
  for(i in 1:length(scramble_res)){
    Ct[,i] = scramble_res[[i]]$count
    Fr[,i] = scramble_res[[i]]$freq
  }
  
  Ct_mean = apply(Ct, 1, mean)
  Fr_mean = apply(Fr, 1, mean)
  
  meanRes = data.frame(
    Source = scramble_res[[1]]$Source,
    Subject = scramble_res[[2]]$Subject, 
    count = Ct_mean, 
    freq = Fr_mean, 
    protocol = scramble_res[[1]]$protocol
  )
  return(meanRes)
}

hmgcr_scramble_rules = lapply(hmgcr_scramble_rules, format_rules)
hmgcr_scramble_res = lapply(hmgcr_scramble_rules, GroupComparator)

hmgcr_scramble_concat = concatScramble(hmgcr_scramble_res)


hscramble_plots = lapply(list(hmgcr_scramble_concat), ResultsPlotter)


hmgcr_scramble_concat$Projection = 'Scramble'
results$hmgcr_giana$Projection = 'HMGCR'

#Aggregated HMGCR Projection vs Scramble Average. 
output = rbind(results$hmgcr_giana, hmgcr_scramble_concat)

spl = split(output, output$Source)

t.test(spl$Biopsy$freq[spl$Biopsy$Projection == 'HMGCR'], 
       spl$Biopsy$freq[spl$Biopsy$Projection == 'Scramble'])

t.test(spl$HC_PBMC$freq[spl$HC_PBMC$Projection == 'HMGCR'], 
       spl$HC_PBMC$freq[spl$HC_PBMC$Projection == 'Scramble'])


t.test(spl$Periphery$freq[spl$Periphery$Projection == 'HMGCR'], 
       spl$Periphery$freq[spl$Periphery$Projection == 'Scramble'])


ggplot(output[output$protocol == 'Adaptive' & output$Source != 'Blood' &  output$Source != 'CMV' ,], 
       aes(x = Source, y = freq * 100, fill = Projection)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(width = 0.3) + 
  ylab('Projection Score') + 
  theme_prism()

hmgcr_agg_output = output

#Aggregated CMV Projection vs Scramble Average. 

# cmv_scramble_rules = lapply(cmv_scramble_rules, format_rules)
# cmv_scramble_res = lapply(cmv_scramble_rules, GroupComparator)
# 
# cmv_scramble_concat = concatScramble(cmv_scramble_res)
# 
# cmv_scramble_concat$Projection = 'Scramble'
# results$cmv_giana$Projection = 'CMV'

#Save all plots. 
for(p in c(1:length(plots))){
  plots[[p]] = lapply(plots[[p]], FUN = function(x){
    return(x + labs(title = names(plots)[p]))
  })
 
  fnames = paste0('./', Sys.Date(), '_', NAME, '_', names(plots)[p], c(1:3))
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

save(hmgcr_agg_output, hmgcr_scramble_rules, rule_sets, file = paste0('./', Sys.Date(), '_ProjectionResults.rda'))