## First let's install and load the ggplot2 and ggfortify libraries:
install.packages('ggplot2')
install.packages('ggfortify')

## Now let's load those libraries so we an ue them:
library('ggplot2')
library('ggfortify')

## Load the TSV from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052:
## Data is already normalizd and log2 transformed:
SCLC_dataset_A <- read.delim(
  "C:/Users/Danie/Desktop/SCLC_project/GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv",
  header= TRUE,
  sep='\t',
  stringsAsFactors = FALSE
)

## Set the rownames to be the column of gene names, then remove
## that column from the dataset
rownames(SCLC_dataset_A) <- SCLC_dataset_A$gene
SCLC_dataset_A$gene <- NULL

## Now let's flip the dataset so the columns represent the genes and the rows
## become the samples. The as.data.frame is to transform the object back into
## a data frame instead of a matrix:
SCLC_dataset_A <- as.data.frame(
  t(SCLC_dataset_A)
)

## Let's do a PCA analysis first:

## The first step is to identify the top 1000 most differential genes.
## This an be done by calculating standard deviation (SD) and identifying the 1000
## genes with the largest SD values:

## Calculate the standard deviation for each column (which represents the SD for each gene
## as the genes are in the columns of the dataset) and set the names of the values
## to be the gene names:
SCLC_SD_values <- apply(
  SCLC_dataset_A,
  2,
  sd
)
names(SCLC_SD_values) <- colnames(SCLC_dataset_A)

## Sort the SD values from largest to smallest:
SCLC_SD_values_sorted <- sort(
  SCLC_SD_values,
  decreasing = TRUE
)

## Identify the top 1000 most differential genes (by highest SD) and get their
## names:
SCLC_top_1000_differential_genes <- SCLC_SD_values_sorted[
  c(1:1000)
]

SCLC_top_1000_differential_gene_names <- names(SCLC_top_1000_differential_genes)

## Create a subset of the dataframe with the top 1000 most differential genes only
SCLC_dataset_A_1000_most_differential <- SCLC_dataset_A[
  ,
  SCLC_top_1000_differential_gene_names
]

## Now let's determine which of the samples are tumor or normal samples
## and create a vector with either 'Normal' or 'Tumor' reflecting whether
## each row in the dataset represents a normal or tumor samples.
## Normal samples have the word 'normal' in their name.
Normal_or_tumor_vector <- ifelse(
  grepl(
    'normal',
    rownames(SCLC_dataset_A)
  ),
  'Normal',
  'Tumor'
)

## Add the normal or tumor vector as a column to the most differential dataset
SCLC_dataset_A_1000_most_differential$Sample_type <- Normal_or_tumor_vector

## Let's create a PCA plot with ggplot2 and ggfortify and plot it:
## Make sure to only plot the 1000 most differential genes and not
## the sample type column:
SCLC_dataset_A_1000_most_differential_pca <- prcomp(
  SCLC_dataset_A_1000_most_differential[
    ,
    c(1:1000)
  ],
  scale. = TRUE
)

suppressWarnings(
  autoplot(
    SCLC_dataset_A_1000_most_differential_pca,
    data= SCLC_dataset_A_1000_most_differential,
    colour= 'Sample_type'
  )
)

## It does look like there are a handful of Tumor samples that seem to cluster
## with the normal samples. Additionally the main axis of variation (PC1) doesn't
## seem to differentiate the bulk of the tumor samples from the normal samples.
## Instead it seems to do more to differentiate some of the tumor samples from
## most of the rest of them.

## Let's identify the genes with the largest difference in expression
## between the tumor and normal samples:

## First let's write a function that reads SCLC_dataset_A and calculates
## a p-value for tumor vs normal expression differential:
T_v_N_t_test_calculator <- function(gene_name){

  ## Get the expression of the gene for the tumor and normal samples:
  tumor_expression <- SCLC_dataset_A[
    !grepl(
      'normal',
      rownames(SCLC_dataset_A)
    ),
    gene_name
  ]

  normal_expression <- SCLC_dataset_A[
    grepl(
      'normal',
      rownames(SCLC_dataset_A)
    ),
    gene_name
  ]

  ## Do a t-test comparing the tumor and normal expression:
  t_test_results <- t.test(
    tumor_expression,
    normal_expression
  )

  ## Return the p-value from the t-test
  return(t_test_results$p.value)

  ## Clear the t-test results:
  rm(t_test_results)
}

## First let's write another function that reads SCLC_dataset_B and calculates
## the mean tumor - mean normal expression difference:
T_v_N_mean_difference_calculator <- function(gene_name){

  ## Get the expression of the gene for the tumor and normal samples:
  tumor_expression <- SCLC_dataset_A[
    !grepl(
      'normal',
      rownames(SCLC_dataset_A)
    ),
    gene_name
  ]

  normal_expression <- SCLC_dataset_A[
    grepl(
      'normal',
      rownames(SCLC_dataset_A)
    ),
    gene_name
  ]

  ##calculate the mean expression of each group:
  tumor_expression_mean <- mean(
    tumor_expression,
    na.rm=TRUE
  )

  normal_expression_mean <- mean(
    normal_expression,
    na.rm = TRUE
  )

  ## Return the difference in the means:
  return(
    (
      tumor_expression_mean - normal_expression_mean
    )
  )
}

## Get the tumor vs. normal p-value for every gene in the dataset
## Using the T_v_N_t_test_calculator function:
T_v_N_p_values <- sapply(
  colnames(SCLC_dataset_A),
  T_v_N_t_test_calculator
)

## Get the mean tumor minus normal expression difference or every gene in the dataset
## Using the T_v_N_t_test_calculator function:
T_v_N_mean_difference_values <- sapply(
  colnames(SCLC_dataset_A),
  T_v_N_mean_difference_calculator
)

## Create a dataframe of information with the t-test p-values, expression differences
## and gene names:
output_df <- data.frame(
  'gene_names'= colnames(SCLC_dataset_A),
  'T_v_N_p_value'= T_v_N_p_values,
  'T_v_N_mean_expression_difference'= T_v_N_mean_difference_values,
  stringsAsFactors = FALSE
)

## Output this dataframe:
write.table(
  output_df,
  'C:/Users/Danie/Desktop/SCLC_project/GSE60052_differential_results_output.tsv',
  sep='\t',
  quote= FALSE
)
