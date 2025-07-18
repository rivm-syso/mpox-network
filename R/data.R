#' Observed age distribution case data from Kamituga
#'
#' This object stores the observed proportion of cases by age group
#' from the Kamituga dataset.
#'
#' @format A data frame with age group and observed proportions.
#'
#' @source Vakaniaki, E.H., Kacita, C., Kinganda-Lusamaki, E. et al.
#' Sustained human outbreak of a new MPXV clade I lineage in eastern
#' Democratic Republic of the Congo. Nat Med 30, 2791–2795 (2024).
#' <https://doi.org/10.1038/s41591-024-03130-3> licensed under a Creative Commons 
#' Attribution 4.0 International License
#'
#' @keywords datasets
"obs_data_Kami"


#' Synthetic contact matrices for community contacts
#'
#' This object stores the synthetic contact matrices for community contacts for 177
#' geographical locations
#'
#' @format A named list of 177 country codes, 16 by 16 numerical matrices.
#'
#' @source Prem K, Zandvoort Kv, Klepac P, Eggo RM, Davies NG, et al. (2024)
#' Projecting contact matrices in 177 geographical regions: An update and comparison
#' with empirical data for the COVID-19 era. PLOS Computational Biology 20(9): e1012454.
#' <https://doi.org/10.1371/journal.pcbi.1012454>
#'
#' @keywords datasets
"contact_all"


#' Synthetic contact matrices for the home location for 177 geographical locations
#'
#' This object stores the synthetic contact matrices for at-home contacts for 177
#' geographical locations
#'
#' @format A named list of 177 country codes, 16 by 16 numerical matrices.
#'
#' @source Prem K, Zandvoort Kv, Klepac P, Eggo RM, Davies NG, et al. (2024)
#' Projecting contact matrices in 177 geographical regions: An update and comparison
#' with empirical data for the COVID-19 era. PLOS Computational Biology 20(9): e1012454.
#' <https://doi.org/10.1371/journal.pcbi.1012454>
#'
#' @keywords datasets
"contact_home"


#' Population age distribution for the DRC
#'
#' Population age distriution for the DRC in age groups of five years with last
#' age group being 100+ years
#'
#' @format tibble with 23 rows and three columns: age_group, male, female
#'
#' @source United Nations, Department of Economic and Social Affairs, Population Division (2024).
#' World Population Prospects 2024, Online Edition. Copyright © 2024 by United Nations,
#' made available under a Creative Commons license CC BY 3.0 IGO
#' <https://population.un.org/wpp/>
#'
#' @keywords datasets
"DRC_pop2024"
