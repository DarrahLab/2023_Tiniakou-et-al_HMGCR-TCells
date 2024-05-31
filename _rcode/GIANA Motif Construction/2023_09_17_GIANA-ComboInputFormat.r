#Format ET Myositis Data and Myositis Periphery TCR Data 
#into input files compatible with GIANA 
#17 September 2023 AAG 

setwd('/media/aag7319/WDRed/ZZZ_RheumAIRR/')
library(dplyr)
library(parallel)
source('./_rfunctions/ProcessAdaptiveFile.R')
load('./data/ET_myositis/combo-hc-myositis.rda')

#Format in accordance with GIANA input requirements. 
#Place cdr3 in first col, TRBV in second col. Ignoring TRBV allele information. 
d = data.frame(cdr3 = alldata$cdr3_aa, 
               trbv = gsub('\\*.*', '', alldata$v_call), 
               sample = alldata$Sample, 
               source = alldata$Source, 
               count = alldata$duplicate_count)
#Filter sequences without V gene or CDR3 amino acid calls. 
d= d[which(d$cdr3 != '' & d$trbv != ''),]

#Expand rows such that each cell/count is one row. 
expanded = data.frame( cdr3 = rep(d$cdr3, d$count),  
                       trbv = rep(d$trbv, d$count),
                       sample = rep(d$sample, d$count), 
                       source = rep(d$source, d$count))
#Filter Myositis data to only include CD154+ T cell populations. 
d1 = expanded[expanded$source == 'Blood',]
d2 = expanded[expanded$source == 'CMV', ]
dh = expanded[expanded$source == 'HC_PBMC',]
d1$trbv = ProcessAdaptiveVgenes(d1$trbv)
d2$trbv = ProcessAdaptiveVgenes(d2$trbv)
dh$trbv = ProcessAdaptiveVgenes(dh$trbv)
#Load Peripheral HMGCR Myositis patient repertoires.------
periphery = mclapply( as.character(c(12211, 13038, 17137, 18124)), mc.cores = 4, FUN = function(id){
  fname = paste0('./data/ET_myositis/2023_09_Periphery/MYO_', id, '_BL.txt')
  d = read.delim(fname, sep = '\t')
  
  d = data.frame(
   cdr3 = d$cdr3aa, 
  #trbv =  gsub('\\*.*', '', d$v), 
   trbv =  gsub('*,.*', '', d$v), #Only include first called v gene. 
  count = d$count
  )
  
  exp = data.frame( cdr3 = rep(d$cdr3, d$count),  
                    trbv = rep(d$trbv, d$count),
                    sample = id, 
                    source = 'Periphery')
  return(exp)
})

d3 = bind_rows(periphery)  #All 4 peripheral blood HMGCR repertoires. 
d3$trbv = ProcessAdaptiveVgenes(d3$trbv)

d4 = bind_rows(d1,d3) #Combined 4 peripheral reps and 5 patients' anti-HMGCR CD154+ TCRs. 
d5 = bind_rows(d2,d3) #Combined 4 peripheral reps and 5 patients' anti-CMV CD154+ TCRs. 
d6 = bind_rows(
  d1[d1$sample %in% c(12211, 13038, 17137, 18124), ], 
  d3
) #Combined 4 peripheral reps and same 4 patients' anti-HMGCR CD154+ TCRs. 

d7 = bind_rows(
  d2[d2$sample %in% c(12211, 13038, 17137, 18124), ], 
  d3
) #Combined 4 peripheral reps and same 4 patients' anti-HMGCR CD154+ TCRs. 

#Saving a combined input of Perpiheral repertoires and either 
#HMGCR/CMV CD154+ populations from all 5 subjects/only the 4 subjects with matched periphery and CD154 stim. 
write.table(d4, 
            file = paste0('./objects/', Sys.Date(), '_GIANA-MyositisInput-HMGCR.txt'), 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE)

write.table(d5, 
            file = paste0('./objects/', Sys.Date(), '_GIANA-MyositisInput-CMV.txt'), 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE)

write.table(d6,
            file = paste0('./objects/', Sys.Date(), '_GIANA-MyositisInputPaired-HMGCR.txt'), 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE)

write.table(d7, 
            file = paste0('./objects/', Sys.Date(), '_GIANA-MyositisInputPaired-CMV.txt'), 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE)

#Lastly, saving RDS with periphery and CD154 populations separately saved. 
saveRDS(d3, file = paste0('./objects/', Sys.Date(), '_HMGCR-periphery.rds'))
saveRDS(d1, file = paste0('./objects/', Sys.Date(), '_HMGCR-CD154HMGCR.rds'))
saveRDS(d2, file = paste0('./objects/', Sys.Date(), '_HMGCR-CD154CMV.rds'))
saveRDS(dh, file = paste0('./objects/', Sys.Date(), '_HC-periphery.rds'))

#CASGLARGSERGSSSTDTQYF