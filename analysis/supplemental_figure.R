############################################################################
# compare the default scenario with the household contact matrix to the
# scenario with the community matrix: time-dependent results
############################################################################
library(mpoxnetwork)
library(patchwork)
library(ggplot2)
library(dplyr)
library(forcats)

# global parameters to the environment
source("analysis/parameters.R")


# default scenario

data_default <- read_simulation_data(
  scenario_name = "baseline",
  data_dir = DATA_DIR
)

df_incidence_default <- prepare_incidence_data_for_plotting(
  data = data_default,
  population_size = population_size,
  age_distribution = age_distribution,
  N_1 = N_1
)

# scenario: community matrix
data_community <- read_simulation_data(
  scenario_name = "community_contact_matrix",
  data_dir = DATA_DIR
)

df_incidence_community <- prepare_incidence_data_for_plotting(
  data = data_community,
  population_size = population_size,
  age_distribution = age_distribution,
  N_1 = N_1
)


###############################################################################
# combine data
df_incidence <- bind_rows(
  list(
    baseline = df_incidence_default,
    community = df_incidence_community
  ),
  .id = "scenario"
) 

data <- bind_rows(
  list(
    baseline = data_default,
    community = data_community
  ),
  .id = "scenario"
) 

plt_cases <- df_incidence |>
  filter(mode == "total") |>
  mutate(age_group = age_group |>
    as.factor() |>
    fct_collapse("0-14" = c(1, 2, 3), "15-49" = c(4, 5, 6, 7, 8, 9, 10), "50+" = 11) |>
    fct_rev()) |>
  group_by(scenario, time, age_group) |>
  summarise(number = sum(number))


############################################################### plot
# incidence by mode of transmission
p1 <- df_incidence |>
  group_by(scenario, mode, time) |>
  summarise(number = sum(number) / population_size) |>
  ggplot(aes(x = time, y = number, col = mode, linetype = scenario)) +
  geom_line(show.legend = c(col = TRUE, linetype = FALSE)) +
  labs(x = "Time (days)", title = "Incidence of infection", y = "", col = "", linetype = "") +
  xlim(c(0, 200)) + 
  theme(legend.position = "bottom")

# time-varying reproduction number
p2 <- data |>
  ggplot(aes(x = time, y = R_t, linetype = scenario)) +
  geom_line() +
  labs(x = "Time (days)", 
       title = "Time-varying\nreproduction number", 
       y = "", 
       linetype = "") +
  xlim(c(0, 200)) +
  ylim(c(0, 1.40)) 

# mean degree of newly infected cases
p3 <- data |>
  ggplot(aes(x = time, y = mean_degree_newly_infected, linetype = scenario)) +
  geom_line() +
  xlim(c(0, 200)) +
  ylim(c(0, 3)) +
  labs(x = "Time (days)", 
       title = "Mean degree of\nnewly infected cases", 
       y = "", 
       linetype = "")

# distribution over age groups: default scenario
p4 <- ggplot(
  data = plt_cases |> filter(scenario == "baseline"),
  aes(x = time, y = number, fill = age_group |> as.factor())
) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  xlim(c(0, 200)) +
  labs(y = "fraction total cases", x = "time", fill = "", title = "Baseline scenario") +
  guides(fill = "none")


# distribution over age groups: community matrix scenario
p5 <- ggplot(
  data = plt_cases |> filter(scenario == "community"),
  aes(x = time, y = number, fill = age_group |> as.factor())
) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  xlim(c(0, 200)) +
  labs(y = "fraction total cases", x = "time", fill = "age group", title = "Community matrix scenario") +
  guides(fill = guide_legend(reverse = TRUE))

# combine plots
p <- (p1 + p2 + p3 +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
) /
  (p4 + p5 +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  plot_annotation(tag_levels = "A")

ggsave(fs::path(DATA_DIR, "supplement_sensitivity_community.png"), 
       plot = p, width = 17.8, height = 12, units = "cm", dpi = 100, bg = "white")

