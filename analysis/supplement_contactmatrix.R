library(mpoxnetwork)
library(patchwork)
library(ggplot2)

source("analysis/parameters.R")


# convert contact matrices for given age groups
# household matrix
contact_matrix_home <- prepare_contact_matrix(contact_home, age_dist, N_1, N_2, M) 
contact_matrix_home <- contact_matrix_home |> as.data.frame()
# community matrix
contact_matrix_all <- prepare_contact_matrix(contact_all, age_dist, N_1, N_2, M)
contact_matrix_all <- contact_matrix_all |> as.data.frame()

# helper vector of age group labels
age_group_labels <- c(
  "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50+"
)
age_group_labels <- factor(age_group_labels, levels = age_group_labels)


# convert matrices to tibbles for plotting purposes
# household matrix
names(contact_matrix_home) <- age_group_labels
contact_matrix_home <- contact_matrix_home |> 
  mutate(index_age_group = age_group_labels,
         type = "household") |> 
  relocate(index_age_group) |> 
  pivot_longer(cols = !c(index_age_group, type), names_to = "contact_age_group", values_to = "contact_rate") |> 
  mutate(contact_age_group = factor(contact_age_group, levels = levels(age_group_labels))) |> 
  arrange(index_age_group, contact_age_group)
  
# community matrix
names(contact_matrix_all) <- age_group_labels
contact_matrix_all <- contact_matrix_all |> 
  mutate(index_age_group = age_group_labels,
         type = "community") |> 
  relocate(index_age_group) |> 
  pivot_longer(cols = !c(index_age_group, type), names_to = "contact_age_group", values_to = "contact_rate") |> 
  mutate(contact_age_group = factor(contact_age_group, levels = levels(age_group_labels))) |> 
  arrange(index_age_group, contact_age_group)

# combine into one tibble
cnt_matrix <- bind_rows(
  contact_matrix_home,
  contact_matrix_all
)

# plot
p_home <- ggplot(
  data = contact_matrix_home,
  mapping = aes(
    x = index_age_group,
    y = contact_age_group,
    size = contact_rate,
    col = type
  )
) +
  coord_fixed() +
  geom_point(col = "#F8766D") +
  scale_size_area() +
  labs(
    x = "index age group",
    y = "contact age group",
    size = "",
    title = "Household"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))

p_communnity <- ggplot(
  data = contact_matrix_all,
  mapping = aes(
    x = index_age_group,
    y = contact_age_group,
    size = contact_rate,
    col = type
  )
) +
  coord_fixed() +
  geom_point(col = "#619CFF") +
  scale_size_area() +
  labs(
    x = "index age group",
    y = "contact age group",
    size = "",
    title = "Community"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))

p <- p_home + p_communnity +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Contact matrices: number of contacts per person per day by age group",
    tag_levels = "A") & 
  theme(legend.position = "bottom")

ggsave(fs::path(DATA_DIR, "contact_matrices.png"), 
       plot = p, width = 22, height = 10, units = "cm", dpi = 100, bg = "white")
