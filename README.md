# HWC

# **Identification of** **cancer** **driver genes based on hierarchical weak consensus model**

Please ensure that R Studio and Matlab are installed before running the code

The following are the code running steps (using BRCA as an example)

1. Open BRCA.rad using R Studio.

2. drivegene.R：Run code from 1 to 120 in R Studio. Export gene_ exp1.txt、mrna_ Exp.csv, means.txt files.

Note：gene_ exp1.txt preserves the threshold and fluctuation values of each gene; mrna_ Exp.csv is the gene expression matrix; means.txt is the mean of each gene.

3. Performed help.m, filtering the matrix to obtain remrna.txt.

Note: The gene expression data of LUAD includes 19429 genes and 572 samples; The gene expression data of PRAD includes 17770 genes and 547 samples.When using these two data files, please pay attention to modifying the size of the matrix.

4. Add the file path of remrna.txt in line 126 of R Studio, and run 126 until the end. Obtaining feature scores

Note：install.packages('igraph')

5. Run fram. m
Note：The cancerene.TXT file is the OMIM benchmark dataset
