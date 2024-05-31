#Heirachical clustering of Eleni Motfis 
#Dropping any motifs which also appear in CMV data. 
#AAG 28 November 2023

setwd('/media/aag7319/WDRed/ZZZ_RheumAIRR/')
library(stringdist)
library(clValid)
#library(dynamicTreeCut)
library(RColorBrewer)
library(circlize)
library(dplyr)
library(dendextend) #color_branches
library(ggseqlogo)
library(msa)
ALPHA = 0.05
hmgcr_cluster_results = read.csv('./results/ET_myositis/2023-09-18_top-giana-clusters.csv')
cmv_cluster_results = read.csv( "./results/ET_myositis/2023-09-18_cmv-cluster-results-filt.csv") 
cmv_rules = data.frame(
  seq = cmv_cluster_results$Consensus[cmv_cluster_results$padj < ALPHA],
  trbv = cmv_cluster_results$vConsensus[cmv_cluster_results$padj < ALPHA])

hmgcr_rules = data.frame(
  seq = hmgcr_cluster_results$Consensus[hmgcr_cluster_results$padj < ALPHA], 
  trbv=  hmgcr_cluster_results$vConsensus[hmgcr_cluster_results$padj < ALPHA])

cmv_rules$uid = paste0(cmv_rules$trbv, '_', cmv_rules$seq)
hmgcr_rules$uid = paste0(hmgcr_rules$trbv, '_', hmgcr_rules$seq)
common = which(hmgcr_rules$uid %in% cmv_rules$uid)

d = hmgcr_cluster_results[hmgcr_cluster_results$padj < ALPHA,][-common,]
excluded = hmgcr_cluster_results[hmgcr_cluster_results$padj < ALPHA,][common,] #these ones are in common w/CMV analysis. 
write.csv(d,row.names = FALSE, file = paste0(Sys.Date(), '_', 'top-giana-clusters-hmgcr-exclusive.csv'))
write.csv(excluded, row.names = FALSE, file = paste0(Sys.Date(), '_', 'top-giana-clusters-hmgcr-common-cmv.csv'))
#Analysis of consensus CDR3 motifs------
motifs = d$Consensus
dists = stringdist::stringdistmatrix(motifs, useNames = 'strings', method = 'lv')

cl = hclust(dists)

dend = as.dendrogram(cl)
# dcl = cutreeDynamic(distM = as.matrix(dists), dendro = cl, minClusterSize = 12, 
#               deepSplit = 1)

# Make Plot Based on Optimal Cluster Number: 8------
cols_set2 = c("#66C2A5" ,"#8DA0CB" ,"#E78AC3", "#A6D854", "#FC8D62" , "#FFD92F",  "#B3B3B3", "#E5C494") 
cols_dark2 = c("#1B9E77" , "#7570B3", "#E7298A", "#66A61E" , "#D95F02", "#E6AB02", "#666666" ,  "#A6761D")
cols = cols_dark2
clusters_final = cutree(cl, k = 8)
dend<- dend %>%
  color_branches(k = 8, col = rev(cols)) %>%
set("branches_lwd", 1.5) 

par(cex=0.3)
plot =  circlize_dendrogram(dend,
                            labels_track_height = 0.3, 
                            dend_track_height = 0.4)

plot
#Generate logo plots for sequences in each heiarchical cluster-----

giana = read.delim('./objects/GIANA_outputs/2023-09-17_GIANA-MyositisInput-HMGCR--RotationEncodingBL62.txt',
                          skip = 2, header = FALSE)

c = cbind( d[d$padj < ALPHA,] , clusters_final)
c = c[order(c$clusters_final, decreasing = FALSE),]
write.csv(c, file = './results/ET_myositis/2023-11-28_top-giana-clusters-byhclust.csv' )

spl = split(c, c$clusters_final)
spl = lapply(spl, FUN = function(x) return(x$cluster)) #GIANA cluster ids for each heiarchical cluster

logos = vector(mode = 'list', length = length(spl))

for(i in 1:length(spl)){
  seqs = giana$V1[giana$V2 %in% spl[[i]]]
  
  x = msa(seqs, method = 'ClustalOmega', type = 'protein')
  x = as.character(x)
  gg = ggplot() + geom_logo(x) + theme_logo()
  logos[[i]] = gg
  
  fname = paste0('./ClusterLogo', i, '.pdf')
  pdf(file = fname, width = 5, height = 2.5)
  plot(gg)
  dev.off()
}

#UMAP to visualize groupings-----
umap = umap(as.matrix(dists))

udf = data.frame(UMAP1 = umap$layout[,1], UMAP2 = umap$layout[,2])

ggplot(udf, aes(x = UMAP1, y = UMAP2)) + 
  geom_point()

#Cluster Performance Metrics using clValid package------
motifs = d$Consensus
dists = stringdist::stringdistmatrix(motifs, useNames = 'strings', method = 'lv')

Report = clValid(as.matrix(dists), nClust = c(2:12))

df = as.data.frame.matrix(t(Report@measures[,,1]))
df$nK = c(2:12)

ggplot(df, aes(x = nK, y = APN )) + geom_path()
ggplot(df, aes(x = nK, y = AD )) + geom_path()
ggplot(df, aes(x = nK, y = ADM )) + geom_path()
ggplot(df, aes(x = nK, y = FOM)) + geom_path()

write.csv(df, file = paste0('./code-outputs/',Sys.Date(), '_cluster-metrics-excluded.csv'))


#Summarize TRBV usage--------
vgenes =  d$vConsensus[d$padj < ALPHA]
table(vgenes) #This is vgene usage in selected clusters. 

hmgcr_result = read.delim('./objects/GIANA_outputs/2023-08-28_GIANA-ComboInput-HMGCR--RotationEncodingBL62.txt',
                          skip = 2, header = FALSE)
colnames(hmgcr_result) = c('seq', 'cluster', 'trbv', 'pid', 'group')
clusters_all = unique(hmgcr_result$V2)
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
system.time(
vgenes_all <- mclapply(clusters_all, dominant_vgene, go = hmgcr_result, mc.cores = 14)

)

vgenes_all = unlist(vgenes_all)

vgenes_allf = table(vgenes_all)