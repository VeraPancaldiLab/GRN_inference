# GRN_inference
Work in progress 

Here we follow the steps from GRN inference of regulatory network (focussing on macrophage polarization) to intracellular Boolean modelling.
The cellular system consist of monocytes and CLL co-culture, which are set to interact for 13 days. From various interactions between the two cell types, monocytes will eventually differentiate and polarize into Nurse Like Cells, which in turn will help the CLL cells to escape apoptosis. 
Our analysis aims to identify the intercellular molecular interactions that establisht the macrophage-CLL interaction in a time-course. For this reason, RNAseq of CLL cells and macrophages are performed in several time-points. Then, we perform extensive network inference and analysis to identify which are the genes and their interactions that explain the transcriptomic dynamics observed in the time-series transcriptomics. To answer this question, we
1. apply inference algorithms for GRN inference. In our case, 2 different methods are applied: dynGENIE3 and SWING-RF. The methods are firstly applied on the CLL cell time-series. 
2. Perform network analysis to indetify the most important genes in the GRN.
3. Perform GO and GSEA to indetify the pathways involved in the CLL cell time-series transcriptomics. 

The same analysis will be perfomed on macrophage time-series transcriptomics. The ultimate goal is to build a mathematical model of macrophage-CLL system, by incorporating both intracellular and intercellular interactions between the 2 cell types, at different states. 
