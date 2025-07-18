### Figure S xx
# load data
result_tibble <- readRDS(file.path(DATA_DIR, "degree_gridsearch_summary.rds"))

filtered_data <- result_tibble %>%
  filter(N_0_val == N_0_val[which.min(abs(N_0_val - N_0))]) %>% # extract the closest value to the default N_0
  filter(k_max < 30) %>%
  filter(2 < alpha) %>%
  filter(alpha < 4)

# set common scale limits
fill_limits <- range(c(filtered_data$R_ss, filtered_data$R_sg), na.rm = TRUE)

# heatmaps
heatmap_Rss_plot <- ggplot(filtered_data, aes(x = alpha_val, y = k_max_val, fill = R_ss)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#005F73", mid = "white", high = "#AE2012", midpoint = 1 # ,
    # limits = fill_limits
  ) +
  labs(
    x = expression("Pareto exponent " * alpha),
    y = expression("Maximum cutoff " * k[max]),
    fill = expression(R["s,s"])
  ) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

heatmap_Rsg_plot <- ggplot(
  filtered_data,
  aes(x = alpha_val, y = k_max_val, fill = R_sg)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#005F73", mid = "white", high = "#AE2012", midpoint = 1 # ,
    # limits = fill_limits
  ) +
  labs(
    x = expression("Pareto exponent " * alpha),
    y = expression("Maximum cutoff " * k[max]),
    fill = expression(R["s,ns"])
  ) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

# full NGM
globalNGM_plot <- ggplot(df_globalNGM, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "#AE2012") +
  labs(x = "Contactor age", y = "Contactee age", fill = "Infection") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

### Figure 1 (Merge all) ----
figure_supp_NGM_plot <- globalNGM_plot | (heatmap_Rss_plot / heatmap_Rsg_plot)
figure_supp_NGM_plot <- figure_supp_NGM_plot +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )

ggsave(fs::path(DATA_DIR, "figure_supp_NGM_plot.png"),
  plot = figure_supp_NGM_plot, width = 17.8, height = 12, units = "cm", dpi = 300, bg = "white"
)
