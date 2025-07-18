# source default parameters
source("analysis/parameters.R")

# shape data frame for plotting NGM submatrices
NGM_result <- NGM_fn(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)

age_groups <- c(
  "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50+"
)
labels <- c(paste0(age_groups, " S"), paste0(age_groups, " G"))

dimnames(NGM_result$k_gg_large) <- list(age_groups, age_groups)
dimnames(NGM_result$k_ss_large) <- list(age_groups, age_groups)
dimnames(NGM_result$t_global_NGM_large) <- list(labels, labels)

df_k_gg <- melt(NGM_result$k_gg_large)
df_k_ss <- melt(NGM_result$k_ss_large)
df_globalNGM <- melt(NGM_result$t_global_NGM_large)

global_NGM_agg <- matrix(
  c(
    max(Re(eigen(NGM_result$k_gs_large)$values)),
    max(Re(eigen(NGM_result$k_gg_large)$values)),
    max(Re(eigen(NGM_result$k_ss_large)$values)),
    max(Re(eigen(NGM_result$k_sg_large)$values))
  ),
  c(2, 2)
)

colnames(global_NGM_agg) <- c("non-sexual", "sexual")
rownames(global_NGM_agg) <- c("sexual", "non-sexual")
df_global_NGM_agg <- melt(global_NGM_agg)
