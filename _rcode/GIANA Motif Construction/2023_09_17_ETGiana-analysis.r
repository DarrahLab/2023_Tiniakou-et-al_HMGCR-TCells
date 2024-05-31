#Re-doing significance tests on GIANA data. 
#Focusing on clusters containing at least on CD154+ TCR from two or more subjects. 
#AAG 28 August 2023

setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
library(dplyr)
library(parallel)
library(umap)
library(ggplot2)
library(mutoss)
NUM='ORIG'
 cmv_result = read.delim('./objects/2023-09-18_GIANA-MyositisInput-CMV--RotationEncodingBL62.txt',
                         skip = 2, header = FALSE)
hmgcr_result = read.delim('./objects/2023-09-17_GIANA-MyositisInput-HMGCR--RotationEncodingBL62.txt',
                          skip = 2, header = FALSE)

# cmv_result = read.delim(paste0('./objects/GIANA_outputs/CMV_', NUM, '.txt'),
#                         skip = 2, header = FALSE)
# hmgcr_result = read.delim(paste0('./objects/GIANA_outputs/HMGCR_', NUM, '.txt'),
#                           skip = 2, header = FALSE)


cmv_input = read.delim('./objects/2023-09-18_GIANA-MyositisInput-CMV.txt',header = TRUE)
hmgcr_input = read.delim('./objects/2023-09-17_GIANA-MyositisInput-HMGCR.txt',header = TRUE)


process_results = function(output, input){
  colnames(output) = c('cdr3_aa', 'cluster', 'trbv', 'subject', 'group')
  total_nseqs = table(input$source)
  contingency_table_template = matrix(data = 0, nrow = 2, ncol = 2)
  contingency_table_template[,1] = as.numeric(total_nseqs)
  rownames(contingency_table_template) = names(total_nseqs)
  
  subjects = unique(output$subject)
  counts =  table(output$cluster, output$subject, output$group) 
  
  
  n_subjects = rowSums(counts[,,1] > 0)
  count_subjects =  rowSums(counts[,,1] )
  n_controls =  rowSums(counts[,,2] > 0)
  count_controls = rowSums(counts[,,2] )
  
  coi = rownames(counts)[n_subjects > 0] #these clusters contain sequences from >=1 HMGCR individual. 
  pvals = unlist(lapply(coi, composition_test, go = output, 
                        contingency_table_template = contingency_table_template, 
                        mode = 'fisher'))
  
  df = data.frame(
    cluster = coi, 
    p = pvals, 
    padj = BY(pvals, alpha = 0.05, silent = TRUE)$adjPValues
  )
  idx = match(df$cluster, rownames(counts))
  df$nseqs = rowSums(counts)[idx]
  df$nSubjects = n_subjects[idx]
  df$nControls = n_controls[idx]
  df$nSeq.Sub =  count_subjects[idx]
  df$nSeq.Ctl =  count_controls[idx]
  df$Consensus = 'Not Computed'
  df$vConsensus =  'Not Computed'
  df$Consensus[df$padj < 0.1] = unlist(lapply(df$cluster[df$padj < 0.1] , process_logo, go = output))
  df$vConsensus[df$padj < 0.1]  = unlist(lapply(df$cluster[df$padj < 0.1] , dominant_vgene, go = output))
  df = df[order(df$padj),]
  rownames(df) = NULL
  return(list(df))
}

process_logo <- function(cluster_id, go, plot = FALSE){
  require(msa)
  require(ggseqlogo)
  
  seqs = go$cdr3_aa[which(go$cluster == cluster_id)]
  
  aligned = msa(seqs, method = 'ClustalOmega', type = 'protein')
  
  con = msaConsensusSequence(aligned)
  
  aligned = as.character(aligned)
  
  if(plot == TRUE){
    gg = ggplot() + geom_logo(aligned) + theme_logo()
    return(list(gg = gg, consensus = con))
  }else{
    return(con)
  }

  
}
dominant_vgene <- function(cluster_id, go){
  vgenes = go$trbv[which(go$cluster == cluster_id)]
  tab = table(vgenes)
  ratio = tab[which(tab == max(tab))] / sum(tab)
  if(length(ratio) > 1){return(NA)}
  
  if(ratio >0.5){
    return(names(tab)[which(tab == max(tab))])
  }else{
    return(NA)
  }
  
}
composition_test <- function(x, go, contingency_table_template, mode = 'fisher'){
  #Takes cluster id as input. 
  #Runs fisher exact to determine statistical significance.  
  cdr3s = go$cdr3_aa[which(go$cluster == x)]
  members = go$group[which(go$cdr3_aa %in% cdr3s)]
  mtab = table(members)
  idx = match(names(mtab), rownames(contingency_table_template))
  
  ctab = contingency_table_template
  ctab[,2] = 0
  ctab[idx,2] = mtab
  ctab[,1] = ctab[,1] - ctab[,2]
  
  if(mode == 'fisher'){
    res = fisher.test(ctab, simulate.p.value = TRUE, workspace = 2e6)
    p = res$p.value
  }else if(mode == 'chisq'){
    res = chisq.test(ctab)
    p = res$p.value
  }
  
  return(p)
}

cmv_clusters = process_results(cmv_result, cmv_input)
hmgcr_clusters = process_results(hmgcr_result,hmgcr_input)

 fname1 = paste0('./code-outputs/', Sys.Date(), '_cmv-cluster-results_', NUM, '.csv')
 write.csv(cmv_clusters, file = fname1, row.names = FALSE)

fname2 = paste0('./code-outputs/', Sys.Date(), '_hmgcr-cluster-results_', NUM, '.csv')
write.csv(hmgcr_clusters, file = fname2, row.names = FALSE)


fname1 = paste0('./code-outputs/', Sys.Date(), '_cmv-cluster-results-filt', '.csv')
write.csv(cmv_cluster_results, file = fname1, row.names = FALSE)

fname2 = paste0('./code-outputs/', Sys.Date(), '_hmgcr-cluster-results-filt_', '.csv')
write.csv(hmgcr_cluster_results, file = fname2, row.names = FALSE)
