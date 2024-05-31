# 2023_Tiniakou-et-al_HMGCR-TCells

This repository contains raw data, processed R objects, and code relevant to the manuscript titled, "Precise Identification and Tracking of HMGCR-reactive CD4+ T Cells in the Target Tissue of Patients with Anti-HMGCR Immune Mediated Necrotizing Myopathy" by Tiniakou et al. Code written and compiled by Alexaxnder Girgis.

For raw TCR repertoire data, please see files in  the folder "/all-tcr-repertoires/" this folder contains: (i) Peripheral TCR repertoire from 4 HMGCR Myositis patients, (ii) Muscle Biopsy repertoires of 5 HMGCR Myositis patients, (iii) Sorted CD4+CD154+ TCRs activated by HMGCR peptides, and (iv) Sorted CD4+CD154+ TCRs activated by CMV peptides. 

Also included are intermediate, processed data which were generated to find TCR motif groups, calculate repertoire characteristics, and perform motif projections. 
When executing code, results files will populate the empty directory /code-outputs/. These files will be identical to intermediate data objects provided in the /objects/ folder. 

Within code, HMGCR+ Myositis subjects are named differently than published aliases used in figures. The key is as follows: 

HMGCR+17==12211
HMGCR+18==13038
HMGCR+19==17137
HMGCR+20==18088
HMGCR+21==18124

Alexander Girgis 
Johns Hopkins University
May 2024

###DESCRIPTION OF CODE FILES#####

####Plot Summaries##### 
2023_05_14_CollectedPlots.rmd				More polished plots from data compiled and saved across previous R scripts. 

####TCR Repertoire Characterization####
2024_05_14_MorisitaOverlap.r				Calculate pairwise Morisita overlap and intersection of repertories. Save output. 
2024_05_14_PeripheryHMGCROverlap.r			Calculature percentage of HMGCR Myositis patient repertoires which contain anti-HMGCR CD4+CD154+ TCRs. Supplemental Figure 3. 
2024_05_14_RDI.r					Calculate RDI between TCR repertoires and save output. 
2024_05_14_ShannonEquitabilityPlot.r			Calculate clonality of TCR repertoires. Supplemental Figure 3. 

####Formation of GIANA Motifs. ####
2023_11_28_GIANAMotifClustering-cmv-omit.r		Find motifs of significant GIANA clusters. Cluster into Heiarchical motif groups. 
2023_09_17_ETGiana-analysis.r				Analyze results of GIANA clustering to identify significant clusters enriched for CD4+ CD451+ sequences. 
2023_09_17_GIANA-ComboInputFormat.r			Create input files for TCR clustering using GIANA v4.1. 
GIANA-execution-script.sh				Shell script used to execute GIANA on the Johns Hopkins Joint High Performance Computing Exchange (JHPCE). 

####Motif Projections####
2023_05_14_KnownTCRProjections.r			Comparison of CD4+CD154+ motifs with VDJDb. 
2024_05_14_GIANA-Result-Projection.r			Projections of HMGCR motifs in aggregate. 
2024_05_14_GIANA-Result-Projection-w6thBiopsy.r		Projections of HMGCR Motifs by heiarchical cluster. 


####Executing Code Locally####
Change the working directory of R scripts to the parent directory of this repository. All code will load relevant data/objects according to relative filepaths as organized in this repository. 
All analyses executed and tested using package versions and hardware listed in SESSIONINFO document, which is the output of print(sessionInfo()) in the R. 
