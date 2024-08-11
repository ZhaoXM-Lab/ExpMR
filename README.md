# An augmented Mendelian randomization approach provides causality of brain imaging features on complex traits in a single biobank-scale dataset
The scripts are for our work which employs a MR framework to infer cell type-specific causal relationships between gene expression and brain-associated complex traits, using eQTL data from eight cell types and large-scale GWASs of 123 imaging-derived phenotypes (IDPs) and 26 brain disorders and behaviors (DBs). 
___
These scripts depend on the following software and R packages:

- R v4.2.0
- GCTA v1.94.0
- SMR v1.3.1
- MetaXcan
- data.table v1.14.4
- tidyr v1.3.0
- PMR v1.0
- doParallel v1.0.17
- foreach v1.5.2
- gsmr v1.1.0
- gridExtra v2.3
- tidyverse v2.0.0
- AnnotationHub v3.6.0
- org.Hs.eg.db v3.18.0
- clusterProfiler v4.10.0
- dplyr v1.1.4
- ggplot2 v3.4.4
- stats v4.3.2
- ggsankey v0.0.9

___
Here, we organized the custom codes of the computational analyses in our study into 3 parts as shown below.

- Part 1. <b>Transcriptome-wide association study (TWAS) analysis</b>.
We conducted TWAS for the GWAS summary statistics of 123 IDPs and 26 DBs to investigate the potential links between gene expression levels and these brain-associated complex traits.

- Part 2. <b>MR_analysis</b>.
We performed two-sample MR analysis for inferring the causal effects of cell type-specific gene expression on brain-associated complex traits using four methods including SMR, Wald radio, PMR-Egger and GSMR. We also performed bidirectional two-sample MR analysis for inferring the relationships between the above DBs and IDPs using PMR-Egger, GSMR, and methods in R package TwoSampleMR.

- Part 3. <b>Follow_up_analysis</b>.
We performed follow-up analysis which included:
(1) merging MR results (1.1.merge MR results.BRVToDB.R, 1.2.merge MR results.CellToBRV.R, 1.3.merge MR results.CellToDB.R, 1.4.merge MR results.CellToWMMP.R, 1.5.merge MR results.DBToBRV.R, 1.6.merge MR results.DBToWMMP.R, 1.7.merge MR results.WMMPToDB.R);
(2) identifying putative causal cell type-specific eQTL target genes (eGenes) for IDPs and DBs (2.1.putative causal eGenes.CellToDB.R, 2.2.putative causal eGenes.CellToIDP.R, 2.3.heatmap of causal eGenes.R).
(3) characterizing the shared causal eGenes among IDPs and DBs (3.1.hypergeometric test for causal eGene overlap.R, 3.2.bar chart and line chart for causal eGene overlap.R, 3.3.upset and venn figure for causal eGene overlap.R, 3.4.pLI score.R, 3.5.gene enrichment analysis.R, 3.6.heatmap of causal eGene ovelap.R).
(4) replication in enternal datasets (4.1.replication rate in BrainMeta.R, 4.2.replication rate in GWAS Catalog.R, 4.3.replication rate in PhychEncode.R).
(5) characterizing the potential causal biological pathways amongst them (5.1.pie chart.R, 5.2.dumbbell chart.BRV.R, 5.3.dumbbell chart.WMMP.R).
(6) exploring their gene expression patterns using external single-cell data (6.1.putative causal DBs.DBToIDP.R, 6.2.putative causal IDPs.IDPToDB.R, 6.3.putative routes.cell_type_eGene-DB-IDP.sankey plot.R).
___
If you have any questions, please contact Anyi Yang (yanganyi_angie@163.com) or Xingzhong Zhao (naturescarl@gmail.com).

