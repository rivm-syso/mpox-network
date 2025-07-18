library(mpoxnetwork)
library(patchwork)
library(ggplot2)
library(dplyr)
library(forcats)

# global parameters to the environment
source("analysis/parameters.R")

df_total <- tibble()

for (scenario_name in c(
  "baseline",
  "beta_0.354703507352208_pcont_0.0971654403777768_sars_0.739422376355021_sarh_0.3",
  "beta_0.250917268462118_pcont_0.151146240587653_sars_0.667480026891617_sarh_0.4",
  "beta_0.125_pcont_0.220833273205528_sars_0.5_sarh_0.493424139165264"
)) {
  data <- read_simulation_data(
    scenario_name = scenario_name,
    data_dir = DATA_DIR
  )

  df_incidence <- prepare_incidence_data_for_plotting(
    data = data,
    population_size = population_size,
    age_distribution = age_distribution,
    N_1 = N_1
  ) |>
    mutate(scenario = scenario_name)

  df_total <- df_total |> bind_rows(
    df_incidence
  )
}

# rename scenarios
df_total <- df_total |>
  mutate(scenario = case_when(
    scenario == "baseline" ~ "baseline",
    scenario == "beta_0.354703507352208_pcont_0.0971654403777768_sars_0.739422376355021_sarh_0.3" ~ "scenario A",
    scenario == "beta_0.250917268462118_pcont_0.151146240587653_sars_0.667480026891617_sarh_0.4" ~ "scenario B",
    scenario == "beta_0.125_pcont_0.220833273205528_sars_0.5_sarh_0.493424139165264" ~ "scenario C"
  )) |>
  mutate(scenario = factor(
    scenario,
    levels = c("baseline", "scenario A", "scenario B", "scenario C")
  ))


############################################################### plot
# incidence by mode of transmission
p1 <- df_total |>
  group_by(scenario, mode, time) |>
  summarise(number = sum(number) / population_size) |>
  ggplot(aes(x = time, y = number, col = mode, linetype = mode, linewidth = mode)) +
  geom_line() +
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
  xlim(c(0, 300)) +
  facet_wrap(~scenario, scales = "free") +
  labs(x = "Time (days)", title = "Incidence of infection", y = "", col = "", linetype = "", linewidth = "") +
  theme(legend.position = "bottom")

ggsave(fs::path(DATA_DIR, "supplement_incidence_simulation_scenarios.png"),
  plot = p1, width = 17.8, height = 11, units = "cm", dpi = 100, bg = "white"
)
