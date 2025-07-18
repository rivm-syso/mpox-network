library(tidyverse)
library(mpoxnetwork)

##############################################################
# final size table for combinations of beta and p_cont
# retrieve data paths
beta_pcont_scenarios <- fs::dir_ls(DATA_DIR)[fs::dir_ls(DATA_DIR) |> stringr::str_detect("data_beta.*.rds$")]

# calculate the final size per mode of transmission and per scenario
tbl_1 <- beta_pcont_scenarios |>
  purrr::map(\(scenario_path) {
    scenario <- fs::path_file(scenario_path)
    scenario <- scenario |> stringr::str_remove(".rds")

    # read appropriate data
    data <- readRDS(scenario_path)
    # incidence of global transmission
    df_incidence_global <- process_data("incidence_global", data, population_size, age_distribution)
    # incidence of sexual transmission
    df_incidence_sexual <- process_data("incidence_sexual", data, population_size, age_distribution, sexual_active_only = TRUE, N_1 = N_1)
    # total incidence over all age groups and modes of transmission
    df_total_incidence <- df_incidence_global |>
      bind_rows(df_incidence_sexual)
    # combine tibbles into one tibble
    df_incidence <- bind_rows(
      list(
        total = df_total_incidence,
        global = df_incidence_global,
        sexual = df_incidence_sexual
      ),
      .id = "mode"
    )
    # plot incidence over time to double check that epidemic has run its course
    p_incidence <- df_incidence |>
      group_by(mode, time) |>
      summarise(number = sum(number) / population_size) |>
      ggplot(aes(x = time, y = number, col = mode, linetype = mode, linewidth = mode)) +
      geom_line()
    ggsave(glue::glue("{DATA_DIR}/plt_incidence_{scenario}.png", file = p_incidence))

    # extract beta and p_cont parameters to add to table
    parms_vals <- scenario_path |> stringr::str_extract_all("(?<=_)[0-9]+\\.[0-9]+")
    beta <- parms_vals[[1]][1] |> as.numeric()
    p_cont <- parms_vals[[1]][2] |> as.numeric()
    sar_s <- parms_vals[[1]][3] |> as.numeric()
    sar_h <- parms_vals[[1]][4] |> as.numeric()

    # (re)-calculate reproduction numbers to add to table
    repro <- reproduction_number(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)
    ngm <- NGM_fn(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)

    # calculate final size
    df_final <- df_incidence |>
      # first sum over age groups by mode over time
      group_by(mode, time) |>
      summarise(number = sum(number)) |>
      # then integrate over the area under the curve for each mode of transmission
      group_by(mode) |>
      # integrate over the area under the curve
      summarise(final = pracma::trapz(time, number)) |>
      # add identifier for each scenario and other helper columns
      mutate(
        beta = beta,
        p_cont = p_cont,
        sar_s = sar_s,
        sar_h = sar_h,
        R_0 = ngm$reproduction_number,
        R_ss = ngm$R_ss,
        R_gg = ngm$R_gg,
        R_sg = ngm$R_sg,
        R_gs = ngm$R_gs,
        final_fraction = final / population_size,
        scenario = scenario
      )
    return(df_final)
  }) |>
  # combine the results into one tibble
  bind_rows()

# pivot table for displaying results: only fractions
tbl_pivot <- tbl_1 |>
  select(-scenario, -final) |>
  tidyr::pivot_wider(names_from = "mode", values_from = c("final_fraction")) |>
  mutate(across(.cols = where(is.numeric), \(x) round(x, 2))) |>
  arrange(beta) |>
  mutate(
    rel = sexual / total,
    rel = rel |> round(2)
  ) 

tbl_pivot <- tbl_pivot |>
  select(-c(beta, p_cont, global, sexual)) |>
  arrange(sar_s |> desc(), sar_h) |> 
  filter(sar_s %in% c(0.8, 0.74, 0.67, 0.5))

tbl_1 |>
  saveRDS(paste0(DATA_DIR, "/final_size_full_table.rds"))


tbl_pivot |>
  saveRDS(paste0(DATA_DIR, "/final_size_results_table.rds"))


tbl_pivot <- readRDS(paste0(DATA_DIR, "/final_size_results_table.rds"))
