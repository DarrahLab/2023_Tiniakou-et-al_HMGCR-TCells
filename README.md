# Precise Identification and Tracking of HMGCR-reactive CD4+ T Cells in the Target Tissue of Patients with Anti-HMGCR Immune Mediated Necrotizing Myopathy

This repository contains raw data, processed R objects, and code relevant to the manuscript titled, "Precise Identification and Tracking of HMGCR-reactive CD4+ T Cells in the Target Tissue of Patients with Anti-HMGCR Immune Mediated Necrotizing Myopathy" by Tiniakou et al. 

Alexander Girgis   
Johns Hopkins University  
May 2024  

## Abstract 

**Background:** Anti-3-hydroxy-3-methylglutaryl coenzyme A reductase (HMGCR)-positive immune mediated necrotizing myopathy (IMNM) is characterized by IgG autoantibodies against HMGCR and a strong association with specific HLA-DR alleles. Although these findings implicate HMGCR-specific CD4+T-cells in disease pathogenesis, no such cells have been described. In this study, we aimed to identify and characterize HMGCR-reactive CD4+T-cells and assess their presence in affected muscle tissue from patients with anti-HMGCR+IMNM.  
  
**Methods:** Peripheral blood mononuclear cells (PBMCs) from patients with anti-HMGCR+IMNM (n=10) and dermatomyositis (DM; n=10) were stimulated with HMGCR protein and peptides identified using a natural antigen processing assay (NAPA; n=6). CD4+T-cell activation was assessed by CD154 upregulation via flow cytometry. T-cell receptor β(TCR) sequencing was performed on paired HMGCR-reactive T-cells and muscle biopsy tissue (n=5).  
  
**Results:** CD4+T-cell responses to HMGCR protein were higher in patients with anti-HMGCR+IMNM compared to DM (median 0.06 vs 0.00, p=0.0059), were enriched in Th1-Th17 cells, and when present, they positively correlated with anti-HMGCR antibody levels (r2=0.89, p=0.0012).  NAPA revealed convergent presentation of seven HMGCR core peptides, with substantial overlap in the peptide repertoires between patients. These HMGCR peptides elicited robust CD4+T-cell responses, with 9/10 anti-HMGCR+IMNM patients responding to at least one peptide, compared to 1/10 DM (p=0.0003). Analysis of HMGCR-reactive TCRs β yielded antigen-reactive motifs that were enriched in muscle biopsies ( projection score 0.03 vs. 0.63, p = 0.007).  
  
**Conclusion:**  HMGCR-antigen-reactive CD4+T-cells are present in the circulation and target tissue of patients with anti-HMGCR+IMNM, correlate with the levels of anti-HMGCR antibodies, and are present in affected muscle tissue, suggesting an active role for these cells in the pathogenesis of anti-HMGCR+IMNM.  


## Raw Repertoire Data 
For raw TCR repertoire data, please see files in  the folder "/all-tcr-repertoires/" this folder contains: (i) Peripheral TCR repertoire from 4 HMGCR Myositis patients, (ii) Muscle Biopsy repertoires of 5 HMGCR Myositis patients, (iii) Sorted CD4+CD154+ TCRs activated by HMGCR peptides, and (iv) Sorted CD4+CD154+ TCRs activated by CMV peptides. 

Also included are intermediate, processed data which were generated to find TCR motif groups, calculate repertoire characteristics, and perform motif projections. 
When executing code, results files will populate the empty directory /code-outputs/. These files will be identical to intermediate data objects provided in the /objects/ folder. 

Within code, HMGCR+ Myositis subjects are named differently than published aliases used in figures. The key is as follows: 
```
HMGCR+17==12211
HMGCR+18==13038
HMGCR+19==17137
HMGCR+20==18088
HMGCR+21==18124
```

## Data Analysis Scripts
To execute code locally, download this repository and change the working directory of R scripts to the parent directory of this repository. All code will load relevant data/objects according to relative filepaths as organized in this repository. All analyses executed and tested using package versions and hardware listed in SESSIONINFO document, which is the output of print(sessionInfo()) in R. 

### Figure Generation
`2023_05_14_CollectedPlots.rmd`				Figures from data compiled and saved across the previous R scripts. 

### TCR Repertoire Characterization
`2024_05_14_MorisitaOverlap.r`				Calculate pairwise Morisita overlap and intersection of repertories. Save output. 
`2024_05_14_PeripheryHMGCROverlap.r`			Calculature percentage of HMGCR Myositis patient repertoires which contain anti-HMGCR CD4+CD154+ TCRs. Supplemental Figure 3.   
`2024_05_14_RDI.r`					Calculate RDI between TCR repertoires and save output.   
`2024_05_14_ShannonEquitabilityPlot.r`			Calculate clonality of TCR repertoires. Supplemental Figure 3.   

### Formation of GIANA Motifs
`2023_11_28_GIANAMotifClustering-cmv-omit.r`		Find motifs of significant GIANA clusters. Cluster into Heiarchical motif groups. 
`2023_09_17_ETGiana-analysis.r`				Analyze results of GIANA clustering to identify significant clusters enriched for CD4+ CD451+ sequences. 
`2023_09_17_GIANA-ComboInputFormat.r`			Create input files for TCR clustering using GIANA v4.1.   
`GIANA-execution-script.sh`				Shell script used to execute GIANA on the Johns Hopkins Joint High Performance Computing Exchange (JHPCE). 

### Motif Projections
`2023_05_14_KnownTCRProjections.r`			Comparison of CD4+CD154+ motifs with VDJDb.   
`2024_05_14_GIANA-Result-Projection.r`			Projections of HMGCR motifs in aggregate.   
`2024_05_14_GIANA-Result-Projection-w6thBiopsy.r`		Projections of HMGCR Motifs by heiarchical cluster.   




