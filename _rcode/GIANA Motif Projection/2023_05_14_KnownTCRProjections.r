#Attempting to feel better about TCR data by cross-checking primary data 
# with clonotypes in VDJDb. I hope to find the following: 
# > CD4+CD154+  anti-CMV T cells have some overlap with known CMV T cells 
# > Not much overlap with CD4+CD154+ anti-HMGCR T cells


#On inspection, all of the above seem to be kinda false. 
#Finding miscellaneous ag reactivity in all data including anti-CMV reaectivity. 
# Also finding SARS-Cov-2 clones which makes no sense as data is all pre-covid. 
#Also finding overlap with vdjdb clones annotated as anti-CD8 despite Eleni's CD154 sorted data 
#which only containing CD4+ T cells. 

#I am therefore placing minimal stock in any of these results, 
# and really the overall quality/annotation rigor of VDJDb period. 
setwd('/media/aag7319/MungoBK/Professional/HMGCR-Tcells/')
library(dplyr)
source('./_rfunctions/ProcessAdaptiveFile.R')
##Load VDJDB----
vdjdb = read.delim('./data/vdjDb-Big-2023-05-10.tsv')
vdjdb$V = ProcessAdaptiveVgenes(vdjdb$V)
vdjdb$V = gsub("\\*.*","",vdjdb$V)
vdjdb$uid = paste0(vdjdb$V, '_', vdjdb$CDR3)
t1 = table(vdjdb$Epitope.species)
t1 = t1[order(t1, decreasing = TRUE)]
print(head(t1, 25)) 

#Load Myositis Data, indluding CD154+ TCRs. ----
load('./objects/combo-hc-myositis.rda')


d1 = data.frame(cdr3 = alldata$cdr3_aa, trbv = alldata$v_call, trbj = alldata$j_call, 
                subject = alldata$Sample, source = alldata$Source, count = alldata$duplicate_count, protocol = 'Adaptive')


combined_data = d1
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
data$uid = paste0(data$trbv, '_', data$cdr3)

#Summarize yields per subject across diff compartments----
blood = as.numeric(table(adaptive_groups$Blood$subject))
cmv = as.numeric(table(adaptive_groups$CMV$subject))
biopsy = as.numeric(table(adaptive_groups$Biopsy$subject))
per = c(as.numeric(table(adaptive_groups$Periphery$subject)[1:3]), 
        NA, as.numeric(table(adaptive_groups$Periphery$subject)[4]))

df = data.frame(
  subject = unique(adaptive_groups$Biopsy$subject), 
 blood = blood, 
 biopsy = biopsy, 
 periphery = per, 
 cmv = cmv
)
write.csv(df,file = './code-outputs/tcr-yields.csv', row.names = FALSE)
#Identify public TCRs among Biopsies -------
b = adaptive_groups$Biopsy
b$trbv = ProcessAdaptiveVgenes(b$trbv)
b$uid = paste0(b$trbv, '_', b$cdr3)
b = b[!is.na(b$cdr3) & b$trbv != 'NA-1*01',]
tab = table(b$uid, b$subject) > 0
ns = rowSums(tab)
t2 = table(ns)

df = data.frame(
  nSubjects = c(1:5), 
  nTCRs = as.numeric(t2)
)
ggplot(df, aes(x = nSubjects, y = nTCRs)) + 
  geom_col() + 
  theme_prism() + 
  scale_y_log10()

 bmat = matrix(data = NA, nrow = 5, ncol = 5) 
 
 spl = split(b, b$subject)
 for( i in 1:5){
   for(j in 1:5){
     bmat[i,j] = sum(spl[[i]]$uid %in% spl[[j]]$uid)
   }
 }

#Clean Workspace----
rm(list = setdiff(ls(), c('data', 'ceft', 'all_hmgcr', 'vdjdb')))
#Helper Functions----
compareVDJDb <- function(input, print_name, db = vdjdb, mode = 'clonotype'){
  if(mode == 'clonotype') matches = which(db$uid %in% input$uid  )
  if(mode == 'cdr3')  matches = which(db$CDR3 %in% input$cdr3 )
  
  if(length(matches) > 0){
    sp = db$Epitope.species[matches]
    print(print_name)
    print(table(sp))
    return(list(db[matches,])) #ouput list length 1 for mapply
  }
  
  return(NA)
}

compareVDJDb2 <- function(input, print_name, db = vdjdb, mode = 'clonotype'){
  #uses grep to perform motif search
  
  if(mode == 'clonotype'){
    arg = paste0(input$uid, collapse = '|')
    matches = grep(arg, db$uid)
  } 
  if(mode == 'cdr3'){
    arg = paste0(input$cdr3, collapse = '|')
    matches = grep(arg, db$CDR3) 
  } 
  
  if(length(matches) > 0){
    sp = db$Epitope.species[matches]
    print(print_name)
    print(table(sp))
    return( list(db[matches,])) #ouput list length 1 for mapply
  }
  
  return(NA)
}
format_rules_motifsearch <- function(rules){
  rules$seq = gsub('\\?', '\\.', rules$seq) #replace unknown chars w/ * for grep search
  rules$seq = gsub('\\%', '\\.', rules$seq)
  rules$seq = gsub('\\*', '\\.', rules$seq)
  
  rules$trbv = gsub('\\*.*', '', rules$trbv) #omit allele specifiers
  
  df = data.frame(
    uid = paste0(rules$trbv, '_', rules$seq), 
    cdr3 = rules$seq
  )
  return(df)
}
#Assess overlap of known CMV/CEFT TCRs with data----
spl = split(data, data$Group)

mode = 'clonotype'
cdr3_overlaps = list(
                     
                      biopsy = compareVDJDb(spl$Biopsy, mode = mode, print_name = 'biopsy'), 
                      cd_hmgcr = compareVDJDb(spl$Blood, mode = mode, print_name = 'cd_hmgcr'), 
                      cd_cmv = compareVDJDb(spl$CMV, mode = mode, print_name = 'cd_cmv'),
                      periphery = compareVDJDb(spl$Periphery, mode = mode, print_name = 'periphery'), 
                      hc.pbmc = compareVDJDb(spl$HC_PBMC, mode = mode, print_name = 'hc.pbmc')
)
                
# mode = 'cdr3'
# cdr3_overlaps = mapply(FUN = compareVDJDb, input =spl, print_name = names(spl), 
#        MoreArgs = list(mode = mode))      
# mode = 'clonotype'
# clone_overlaps =  mapply(FUN = compareVDJDb, input =spl, print_name = names(spl), 
#                                         MoreArgs = list(mode = mode))   
#Assess overlap of CMV/HMGCR Motifs, rather than indiv TCRs, with data-----
#First formatting rules for grep search. 
hmgcr_cluster_results = read.csv("./results/ET_myositis/2023-09-18_hmgcr-cluster-results-filt_.csv") #GIANA
cmv_cluster_results = read.csv( "./results/ET_myositis/2023-09-18_cmv-cluster-results-filt.csv") 
load("./objects/ET-tcr/2023-09-18_Scrambled-ClusterResults.rda") #Null distributions


hmgcr_rules = data.frame(seq = hmgcr_cluster_results$Consensus[hmgcr_cluster_results$padj < ALPHA],
                         trbv = hmgcr_cluster_results$vConsensus[hmgcr_cluster_results$padj < ALPHA])     

idx = match(hmgcr_rules$seq, c$Consensus)
hmgcr_rules_by_clust = split(hmgcr_rules, c$clusters_final[idx])
names(hmgcr_rules_by_clust) = paste0('Group', c(1:8))


cmv_rules = list(cmv = data.frame(
  seq = cmv_cluster_results$Consensus[cmv_cluster_results$padj < ALPHA],
  trbv = cmv_cluster_results$vConsensus[cmv_cluster_results$padj < ALPHA]))

rule_sets = c(hmgcr_rules_by_clust, cmv_rules)

rule_sets = lapply(rule_sets, format_rules_motifsearch)
#Now doing the searches 

# mode = 'cdr3'
# cdr3_overlaps = mapply(FUN = compareVDJDb2, input =rule_sets, print_name = names(rule_sets), 
#                        MoreArgs = list(mode = mode))
mode = 'clonotype'
clone_overlaps = mapply(FUN = compareVDJDb2, input =rule_sets, print_name = names(rule_sets), 
                        MoreArgs = list(mode = mode))
