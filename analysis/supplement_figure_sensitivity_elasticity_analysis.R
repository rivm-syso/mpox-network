library(patchwork)
##############################################################
# compare elasticity and sensitivity of mortality-weighted NGM
# elasticity and sensitivity of NGM at t = 0 and 200
result_pertub_snapshots_long_alt <- readRDS(fs::path(DATA_DIR, "sensitivity_result_pertub_snapshots_long.rds"))
result_pertub_snapshots_long <- readRDS(fs::path(DATA_DIR, "result_pertub_snapshots_long.rds"))

result <- dplyr::bind_rows(
  list(
    baseline = result_pertub_snapshots_long,
    alternative = result_pertub_snapshots_long_alt
  ), .id = "scenario")

result <- result |> dplyr::mutate(scenario = factor(scenario, levels = c("baseline", "alternative")))


ela_plot_0 <- ggplot(
  result |>
    filter(time == 0, type == "elasticity", label == "Mortality"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = scenario)
) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0), labels = c(0.2, 0.1, 0)) + # Show values as positive
  scale_fill_manual(
    values = c("baseline" = "#005F73", "alternative" = "#b2abd2"),
    labels = c("baseline" = "Baseline", "alternative"= "Alternative")
  ) +
  labs(
    x = "Elasticity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  labs(tag = "A")

sens_plot_0 <- ggplot(
  result |>
    filter(time == 0, type == "elasticity", label == "Mortality"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = scenario)
) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  scale_x_continuous(labels = abs) + # Show values as positive
  scale_fill_manual(
    values = c("baseline" = "#94D2BD", "alternative" = "#e66101"),
    labels = c("baseline" = "Baseline", "alternative"= "Alternative")
  ) +
  labs(
    x = "Sensitivity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  labs(tag = "B")

ela_plot_200 <- ggplot(
  result |>
    filter(time == 200, type == "elasticity", label == "Mortality"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = scenario)
) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0), labels = c(0.2, 0.1, 0)) + # Show values as positive
  scale_fill_manual(
    values = c("baseline" = "#005F73", "alternative" = "#b2abd2"),
    labels = c("baseline" = "Baseline", "alternative"= "Alternative")
  ) +
  labs(
    x = "Elasticity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  labs(tag = "C")

sens_plot_200 <- ggplot(
  result |>
    filter(time == 200, type == "elasticity", label == "Mortality"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = scenario)
) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  scale_x_continuous(labels = abs) + # Show values as positive
  scale_fill_manual(
    values = c("baseline" = "#94D2BD", "alternative" = "#e66101"),
    labels = c("baseline" = "Baseline", "alternative"= "Alternative")
  ) +
  labs(
    x = "Sensitivity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  labs(tag = "D")

p1 <- ela_plot_0 + sens_plot_0 + plot_annotation(title = "Early stages of the epidemic") 
p2 <- ela_plot_200 + sens_plot_200 + plot_annotation(title = "Later stage of the epidemic")
figure <- wrap_elements(p1) / wrap_elements(p2) 

ggsave(fs::path(DATA_DIR, "figure_supplement_alternative_mortality.png"),
       plot = figure, width = 13, height = 15, units = "cm", dpi = 100, bg = "white"
)

