library(tidyverse)
library(patchwork)
library(mpoxnetwork)


# load parameters into non-sexual environment
source("analysis/parameters.R")


data <- read_simulation_data(
  scenario_name = "baseline",
  data_dir = DATA_DIR
)

df_incidence <- prepare_incidence_data_for_plotting(
  data = data,
  population_size = population_size,
  age_distribution = age_distribution,
  N_1 = N_1
)

result_pertub_snapshots_long <- readRDS(fs::path(DATA_DIR, "result_pertub_snapshots_long.rds"))

###############################
# restrict times to where there is at least one age group with at least one case
plt_cases <- df_incidence |>
  filter(mode == "total") |>
  mutate(age_group = age_group |>
    as.factor() |>
    fct_collapse("0-14" = c(1, 2, 3), "15-49" = c(4, 5, 6, 7, 8, 9, 10), "50+" = 11) |>
    fct_rev()) |>
  group_by(time, age_group) |>
  summarise(number = sum(number))

time_restricted_cases <- plt_cases |>
  group_by(time) |>
  filter(number > 1) |>
  distinct(time) |>
  pull(time)
###################################

# plot incidence summed by age groups over time stratified by mode of transmission and overall
p_incidence <- df_incidence |>
  group_by(mode, time) |>
  summarise(number = sum(number) / population_size) |>
  ggplot(aes(x = time, y = number, col = mode, linetype = mode, linewidth = mode)) +
  geom_line() +
  geom_line(
    data = data,
    aes(
      x = time,
      y = R_t / 800,
      col = "time-varying reproduction number",
      linetype = "time-varying reproduction number",
      linewidth = "time-varying reproduction number"
    )
  ) +
  scale_color_manual(
    values = c("non_sexual" = "black", "sexual" = "black", "total" = "darkgreen", "time-varying reproduction number" = "orange"),
    breaks = c("non_sexual", "sexual", "total", "time-varying reproduction number"),
    labels = c("non-sexual", "sexual", "total", "time-varying reproduction number")
  ) +
  scale_linetype_manual(
    values = c("non_sexual" = 3, "sexual" = 1, "total" = 1, "time-varying reproduction number" = 4),
    breaks = c("non_sexual", "sexual", "total", "time-varying reproduction number"),
    labels = c("non-sexual", "sexual", "total", "time-varying reproduction number")
  ) +
  scale_linewidth_manual(
    values = c("non_sexual" = 1, "sexual" = 1, "total" = 1.5, "time-varying reproduction number" = 1.5),
    breaks = c("non_sexual", "sexual", "total", "time-varying reproduction number"),
    labels = c("non-sexual", "sexual", "total", "time-varying reproduction number")
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 800, name = "Time-varying\n reproduction number")) +
  labs(x = "Time (days)", y = "Incidence of infection", col = "", linewidth = "", linetype = "") +
  xlim(c(0, 200)) +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))

p_incidence

################################################################################
# percent stacked barchart over time across age groups
# mean degree of newly infected individuals over time


plt_cases <- plt_cases |>
  # put number of cases as NA when the total number of cases is not at least one
  mutate(number = if_else(time %in% time_restricted_cases, number, NA))

# join with mean degree for plotting and put the same restriction on plotting times
plt_cases <- plt_cases |>
  left_join(data |> select(time, mean_degree_newly_infected),
    by = join_by(time)
  ) |>
  mutate(mean_degree_newly_infected = if_else(time %in% time_restricted_cases, mean_degree_newly_infected, NA))

p_distribution_age <- ggplot(
  data = plt_cases,
  aes(x = time, y = number, fill = age_group |> as.factor())
) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  geom_line(aes(x = time, y = mean_degree_newly_infected / 4), col = "red", linetype = 3, linewidth = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 4, name = "mean degree", breaks = 0:10)) +
  labs(y = "fraction total cases", x = "time", fill = "age group") +
  xlim(c(0, 200)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(reverse = TRUE))
  
p_distribution_age

##############################################################
# elasticity and sensitivity of NGM at t = 0 and 200

ela_plot_0 <- ggplot(
  result_pertub_snapshots_long |>
    filter(time == 0, type == "elasticity"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = label)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0), labels = c(0.2, 0.1, 0)) + # Show values as positive
  scale_fill_manual(
    values = c("Infection" = "#BB3E03", "Mortality" = "#005F73"),
    labels = c("Infection", "Mortality")
  ) +
  labs(
    x = "Elasticity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

sens_plot_0 <- ggplot(
  result_pertub_snapshots_long |>
    filter(time == 0, type == "sensitivity"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = label)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(labels = abs) + # Show values as positive
  scale_fill_manual(
    values = c("Infection" = "#EE9B00", "Mortality" = "#94D2BD"),
    labels = c("Infection", "Mortality")
  ) +
  labs(
    x = "Sensitivity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ela_plot_200 <- ggplot(
  result_pertub_snapshots_long |>
    filter(time == 200, type == "elasticity"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = label)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0), labels = c(0.2, 0.1, 0)) + # Show values as positive
  scale_fill_manual(
    values = c("Infection" = "#BB3E03", "Mortality" = "#005F73"),
    labels = c("Infection", "Mortality")
  ) +
  labs(
    x = "Elasticity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

sens_plot_200 <- ggplot(
  result_pertub_snapshots_long |>
    filter(time == 200, type == "sensitivity"),
  aes(x = value, y = factor(group, levels = unique(group)), fill = label)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(labels = abs) + # Show values as positive
  scale_fill_manual(
    values = c("Infection" = "#EE9B00", "Mortality" = "#94D2BD"),
    labels = c("Infection", "Mortality")
  ) +
  labs(
    x = "Sensitivity",
    y = "Age group",
    fill = ""
  ) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#######################################################################
panel_AB <- p_incidence / p_distribution_age

panel_CD <- (ela_plot_0 | sens_plot_0)
panel_EF <- (ela_plot_200 | sens_plot_200)

figure <- panel_AB | (panel_CD / panel_EF) 

figure <- figure + plot_annotation(tag_levels = "A")

ggsave(fs::path(DATA_DIR, "figure_2.png"),
  plot = figure, width = 17.8, height = 15, units = "cm", dpi = 100, bg = "white"
)
