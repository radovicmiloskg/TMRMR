TMRMR package Version 1.1 (June. 17, 2016)

TMRMR: Minimum redundancy maximum relevance feature selection approach for temporal gene expression data.

File geneSelection_KFold.m can be use to reproduce results from the paper. Results are written in the folder RESULTS.

File TMRMR.m can be used to perform TMRMR-C and TMRMR-M feature selection for 'H3N2', 'Rhino' or 'RSV' datasets. Results are written in the folder RESULTS.

File GeneSelection.m can be used to perform feature selection by mRMR, F-statistic, RELIEFF, MT-LASSO, TMRMR-C and TMRMR-M algorithms for 'H3N2', 'Rhino' or 'RSV' datasets. Results are written in the folder RESULTS.

File GeneSelectionNewData.m can be used to perform feature selection by mRMR, F-statistic, RELIEFF, MT-LASSO, TMRMR-C and TMRMR-M algorithms for datasets other than ones used in the paper. New dataset should be present within the working directory as Dataset.mat file. This file should contain the following variables:
1) DATA - GxTxN matrix
2) LABELS - Nx1 vector
where N is the number of subjects, G is the number of genes and T is the number of timesteps.
