## First let's install and load the ggplot2 and ggfortify libraries:
install.packages('ggplot2')
install.packages('ggfortify')

## Now let's load those libraries so we an ue them:
library('ggplot2')
library('ggfortify')

## Load the microarray data from 18 tumor and 18 normal samples from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149507
SCLC_dataset_B <- read.delim(
  "C:/Users/Danie/Desktop/SCLC_project/GSE149507_geo_expr.csv.gz",
  sep=',',
  stringsAsFactors = FALSE
)

## Set the rownames to be the gene names:
## The gene_names are not unique, so try pasting together the EntrezID
## and the gene names to make them unique:
rownames(SCLC_dataset_B) <- paste(
  SCLC_dataset_B$Symbol,
  SCLC_dataset_B$EntrezID,
  sep='_'
)

## Now remove the EntrezID and the gene name from the dataframe:
SCLC_dataset_B$EntrezID <- NULL
SCLC_dataset_B$Symbol <- NULL

## Let's transform the dataset reversing the rows and columns
## (So genes in the columns and samples in the rows):
SCLC_dataset_B <- as.data.frame(
  t(
    SCLC_dataset_B
  )
)

## Calculate the standard deviation for each column (which represents the SD for each gene
## as the genes are in the columns of the dataset) and set the names of the values
## to be the gene names:
SCLC_SD_values_B <- apply(
  SCLC_dataset_B,
  2,
  sd
)
names(SCLC_SD_values_B) <- colnames(SCLC_dataset_B)

## Sort the SD values from largest to smallest:
SCLC_SD_values_sorted_B <- sort(
  SCLC_SD_values_B,
  decreasing = TRUE
)

## Identify the top 1000 most differential genes (by highest SD) and get their
## names:
SCLC_top_1000_differential_genes_B <- SCLC_SD_values_sorted_B[
  c(1:1000)
]

SCLC_top_1000_differential_gene_names_B <- names(SCLC_top_1000_differential_genes_B)

## Create a subset of the dataframe with the top 1000 most differential genes only
SCLC_dataset_B_1000_most_differential <- SCLC_dataset_B[
  ,
  SCLC_top_1000_differential_gene_names_B
]

## Now let's determine which of the samples are tumor or normal samples
## and create a vector with either 'Normal' or 'Tumor' reflecting whether
## each row in the dataset represents a normal or tumor samples.
## Normal samples have the string '_n' in their name.
Normal_or_tumor_vector_B <- ifelse(
  grepl(
    '_n',
    rownames(SCLC_dataset_B)
  ),
  'Normal',
  'Tumor'
)

## Add the normal or tumor vector as a column to the most differential dataset
SCLC_dataset_B_1000_most_differential$Sample_type <- Normal_or_tumor_vector_B

## Let's create a PCA plot with ggplot2 and ggfortify and plot it:
## Make sure to only plot the 1000 most differential genes and not
## the sample type column:
SCLC_dataset_B_1000_most_differential_pca <- prcomp(
  SCLC_dataset_B_1000_most_differential[
    ,
    c(1:1000)
  ],
  scale. = TRUE
)

suppressWarnings(
  autoplot(
    SCLC_dataset_B_1000_most_differential_pca,
    data= SCLC_dataset_B_1000_most_differential,
    colour= 'Sample_type'
  )
)

## Let's identify the genes with the largest difference in expression
## between the tumor and normal samples:

## First let's write a function that reads SCLC_dataset_B and calculates
## a p-value for tumor vs normal expression differential:
T_v_N_t_test_calculator <- function(gene_name){

  ## Get the expression of the gene for the tumor and normal samples:
  tumor_expression <- SCLC_dataset_B[
    grepl(
      '_ca',
      rownames(SCLC_dataset_B)
    ),
    gene_name
  ]

  normal_expression <- SCLC_dataset_B[
    grepl(
      '_n',
      rownames(SCLC_dataset_B)
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
  tumor_expression <- SCLC_dataset_B[
    grepl(
      '_ca',
      rownames(SCLC_dataset_B)
    ),
    gene_name
  ]

  normal_expression <- SCLC_dataset_B[
    grepl(
      '_n',
      rownames(SCLC_dataset_B)
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
  colnames(SCLC_dataset_B),
  T_v_N_t_test_calculator
)

## Get the mean tumor minus normal expression difference or every gene in the dataset
## Using the T_v_N_t_test_calculator function:
T_v_N_mean_difference_values <- sapply(
  colnames(SCLC_dataset_B),
  T_v_N_mean_difference_calculator
)

## Create a dataframe of information with the t-test p-values, expression differences
## and gene names:
output_df <- data.frame(
  'unique_gene_IDs'= colnames(SCLC_dataset_B),
  'gene_names'= sub('\\_.*', '', colnames(SCLC_dataset_B)),
  'T_v_N_p_value'= T_v_N_p_values,
  'T_v_N_mean_expression_difference'= T_v_N_mean_difference_values,
  stringsAsFactors = FALSE
)

## Output this dataframe:
write.table(
  output_df,
  'C:/Users/Danie/Desktop/SCLC_project/GSE149507_differential_results_output.tsv',
  sep='\t',
  quote= FALSE
)
