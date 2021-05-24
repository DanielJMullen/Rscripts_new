## Script to identify lncRNA genes of interest that are present in Evelyn Tran's
## RNA-seq data:

## Data originally sent in email 5/24/2021 from Michele Ramos Correa

## Load devtools
if (!require("devtools", quietly = TRUE)){

  install.packages("devtools")

}

## Load TENETR.data
if (!require("TENETR.data", quietly = TRUE)){

  devtools::install_github("DanielJMullen/TENETR.data")

}

gtf <- TENETR.data::gencode_v22_annotations

gencode_v22_genes <- gtf[
  gtf$type=='gene',
]

rownames(gencode_v22_genes) <- sub(
  '\\..*',
  '',
  gencode_v22_genes$gene_id
)

## Read the first file into R:
## Note: although the file is given as a .csv, this appears to actually be a
## .tsv file:
Evelyn_gene_expression <- read.table(
  "C:/Users/Danie/Downloads/AEC_lines-HLF-Primary_AECs_raw_counts_table.csv",
  sep='\t',
  stringsAsFactors = FALSE,
  header= TRUE
)

## Set the row names in Evelyn's RNAseq data to be the gene names:
rownames(Evelyn_gene_expression) <- Evelyn_gene_expression$Gene_Symbol

lnc_RNAs <- read.table(
  "C:/Users/Danie/Downloads/93_lncRNAs.csv",
  sep=',',
  stringsAsFactors = FALSE,
  header= TRUE
)

## Get gene names for the ENSGs in the lnc RNA data:
lnc_RNA_data_gene_names <- gencode_v22_genes[
  lnc_RNAs$X,
  'gene_name'
]

## Now that we have gene names, check the data in Evelyn's RNA-seq data:
Evelyn_gene_expression_lnc_RNAs <- Evelyn_gene_expression[
  lnc_RNA_data_gene_names,
]

## Lots of genes weren't found, so let's reduce the dataset down to the ones that are:
Evelyn_gene_expression_lnc_RNAs_present <- Evelyn_gene_expression_lnc_RNAs[
  complete.cases(Evelyn_gene_expression_lnc_RNAs),
]

