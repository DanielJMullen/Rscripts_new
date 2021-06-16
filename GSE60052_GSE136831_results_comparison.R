GSE60052_results <- read.delim(
  'C:/Users/Danie/Desktop/SCLC_project/GSE60052_differential_results_output.tsv',
  header= TRUE,
  sep='\t',
  stringsAsFactors = FALSE
)

GSE149507_results <- read.delim(
  'C:/Users/Danie/Desktop/SCLC_project/GSE149507_differential_results_output.tsv',
  header= TRUE,
  sep='\t',
  stringsAsFactors = FALSE
)

## The two datasets analyzed different genes, so let's identify the genes
## present in both:
genes_in_both <- intersect(
  GSE60052_results$gene_names,
  GSE149507_results$gene_names
)

## Let's subset the results of both datasets to include only the genes found
## in each:
GSE60052_results_both <- GSE60052_results[
  genes_in_both,
]

GSE149507_results_both <- GSE149507_results[
  GSE149507_results$gene_names%in%genes_in_both,
]

## Do an FDR correction on the genes found in both datasets:
GSE60052_results_both$T_v_N_p_value_FDR <- p.adjust(
  GSE60052_results_both$T_v_N_p_value,
  method='fdr'
)

GSE149507_results_both$T_v_N_p_value_FDR <- p.adjust(
  GSE149507_results_both$T_v_N_p_value,
  method='fdr'
)

## Identify the counts of genes with FDR p-values <0.05 and FC >1:
nrow(
  GSE60052_results_both[
    GSE60052_results_both$T_v_N_p_value_FDR<0.05,
  ]
)

nrow(
  GSE60052_results_both[
    GSE60052_results_both$T_v_N_mean_expression_difference>1,
  ]
)

nrow(
  GSE149507_results_both[
    GSE149507_results_both$T_v_N_p_value_FDR<0.05,
  ]
)

nrow(
  GSE149507_results_both[
    GSE149507_results_both$T_v_N_mean_expression_difference>1,
  ]
)


## For each gene in each dataset lets calculate a rank both for most significant
## p-value and largest tumor - nomral difference (positive):

## Sort the two datasets by most significant (smallest) p-value:
GSE60052_results_both_p_value_sorted <- GSE60052_results_both[
  order(
    GSE60052_results_both$T_v_N_p_value,
    decreasing= FALSE
  ),
]

GSE149507_results_both_p_value_sorted <- GSE149507_results_both[
  order(
    GSE149507_results_both$T_v_N_p_value,
    decreasing= FALSE
  ),
]

## Sort the two datasets by largest (positive) tumor minus normal mean
## expression difference
GSE60052_results_both_expression_diff_sorted <- GSE60052_results_both[
  order(
    GSE60052_results_both$T_v_N_mean_expression_difference,
    decreasing= TRUE
  ),
]

GSE149507_results_both_expression_diff_sorted <- GSE149507_results_both[
  order(
    GSE149507_results_both$T_v_N_mean_expression_difference,
    decreasing= TRUE
  ),
]

## Change rownames for the GSE149507 datasets to be the gene names only:
rownames(GSE149507_results_both_p_value_sorted) <- GSE149507_results_both_p_value_sorted$gene_names
rownames(GSE149507_results_both_expression_diff_sorted) <- GSE149507_results_both_expression_diff_sorted$gene_names

## Now to each of the sorted datasets, add a rank to each gene from 1 to the
## number of genes in the dataset:
GSE60052_results_both_p_value_sorted$rank <- c(
  1:nrow(GSE60052_results_both_p_value_sorted)
)

GSE149507_results_both_p_value_sorted$rank <- c(
  1:nrow(GSE149507_results_both_p_value_sorted)
)

GSE60052_results_both_expression_diff_sorted$rank <- c(
  1:nrow(GSE60052_results_both_expression_diff_sorted)
)

GSE149507_results_both_expression_diff_sorted$rank <- c(
  1:nrow(GSE149507_results_both_expression_diff_sorted)
)

## Now for each gene let's record the p-value, differential mean expression,
## and rank for each values in the two datasets:
GSE60052_GSE149507_results_df <- data.frame(
  'GSE60052_p_value'= GSE60052_results_both_p_value_sorted[
    genes_in_both,
    'T_v_N_p_value'
  ],
  'GSE60052_p_value_FDR'= GSE60052_results_both_p_value_sorted[
    genes_in_both,
    'T_v_N_p_value_FDR'
  ],
  'GSE60052_p_value_rank'= GSE60052_results_both_p_value_sorted[
    genes_in_both,
    'rank'
  ],
  'GSE60052_TvN_expression_differential'= GSE60052_results_both_expression_diff_sorted[
    genes_in_both,
    'T_v_N_mean_expression_difference'
  ],
  'GSE60052_TvN_expression_differential_rank'= GSE60052_results_both_expression_diff_sorted[
    genes_in_both,
    'rank'
  ],
  'GSE149507_p_value'= GSE149507_results_both_p_value_sorted[
    genes_in_both,
    'T_v_N_p_value'
  ],
  'GSE149507_p_value_FDR'= GSE149507_results_both_p_value_sorted[
    genes_in_both,
    'T_v_N_p_value_FDR'
  ],
  'GSE149507_p_value_rank'= GSE149507_results_both_p_value_sorted[
    genes_in_both,
    'rank'
  ],
  'GSE149507_TvN_expression_differential'= GSE149507_results_both_expression_diff_sorted[
    genes_in_both,
    'T_v_N_mean_expression_difference'
  ],
  'GSE149507_TvN_expression_differential_rank'= GSE149507_results_both_expression_diff_sorted[
    genes_in_both,
    'rank'
  ],
  stringsAsFactors = FALSE
)

## Set the rownames:
rownames(GSE60052_GSE149507_results_df) <- genes_in_both

## Now that we have ranks from both, let's sum the two types of ranks from the
## datasets and see which have the lowest composite ranks:
GSE60052_GSE149507_results_df$p_value_ranks_sum <- GSE60052_GSE149507_results_df$GSE60052_p_value_rank+GSE60052_GSE149507_results_df$GSE149507_p_value_rank
GSE60052_GSE149507_results_df$TvN_expression_differential_ranks_sum <- GSE60052_GSE149507_results_df$GSE60052_TvN_expression_differential_rank+GSE60052_GSE149507_results_df$GSE149507_TvN_expression_differential_rank

## Check the number of genes FDR significant in each:
nrow(
  GSE60052_GSE149507_results_df[
    (
      (GSE60052_GSE149507_results_df$GSE60052_p_value_FDR<0.05) & (GSE60052_GSE149507_results_df$GSE149507_p_value_FDR<0.05)
    ),
  ]
)

nrow(
  GSE60052_GSE149507_results_df[
    (
      (GSE60052_GSE149507_results_df$GSE60052_TvN_expression_differential>1) & (GSE60052_GSE149507_results_df$GSE149507_TvN_expression_differential>1)
    ),
  ]
)

## Get the genes with the lowest rank but with positive FC:
top_genes_by_p_value_rank_positive <- GSE60052_GSE149507_results_df[
  GSE60052_GSE149507_results_df$GSE60052_TvN_expression_differential>0,
]

top_genes_by_p_value_rank_positive <- top_genes_by_p_value_rank_positive[
  order(
    top_genes_by_p_value_rank_positive$p_value_ranks_sum,
    decreasing = FALSE
  ),
]

## Get the genes by top positive FC:
top_genes_by_TvN_expression_differential_rank_positive <- GSE60052_GSE149507_results_df[
  GSE60052_GSE149507_results_df$GSE149507_TvN_expression_differential>0,
]

top_genes_by_TvN_expression_differential_rank_positive <- top_genes_by_TvN_expression_differential_rank_positive[
  order(
    top_genes_by_TvN_expression_differential_rank_positive$TvN_expression_differential_ranks_sum,
    decreasing = FALSE
  ),
]
