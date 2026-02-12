# Environment setup
base::rm(list = base::ls()) 
base::library(tidyverse)
base::library(readxl)
base::library(ggplot2)
base::library(maps)
base::library(scales)

# Import data
data_path <- "../2_Scripts/clean_self_ref_paper_excluded_age.csv"

# ====================================================================
# Figure S1: Geographical visulization 
# ====================================================================

# Part 1: Map Data Preparation
# Fetch world map data using ggplot2 and filter with dplyr
world <- ggplot2::map_data("world") %>% 
  dplyr::filter(region != "Antarctica")

# Import sample metadata using readr
Sample_Info <- readr::read_csv(data_path, col_names = TRUE)

# Aggregate sample sizes by country using dplyr
country_summary <- Sample_Info %>% 
  dplyr::group_by(Data_Collection_country) %>% 
  dplyr::summarise(Total_Sample = base::sum(Sample_patients, na.rm = TRUE))

# Merge map coordinates with statistical data using dplyr
world_with_data <- world %>% 
  dplyr::left_join(country_summary, by = c("region" = "Data_Collection_country"))

# Generate World Map Plot
p_map <- ggplot2::ggplot() +
  # Draw polygons and map fill color to sample size
  ggplot2::geom_polygon(
    data = world_with_data,
    ggplot2::aes(x = long, y = lat, group = group, fill = Total_Sample),
    color = "gray50", linewidth = 0.3
  ) +
  # Define color gradient for the fill
  ggplot2::scale_fill_gradient(
    low = "#cfe7c4",  # Light green for low counts
    high = "#4f845c", # Dark green for high counts
    na.value = "white",
    name = "Sample Size"
  ) +
  # Set map projection and clean theme
  ggplot2::coord_quickmap() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = "bottom"
  ) +
  ggplot2::labs(title = NULL, x = NULL, y = NULL)

base::print(p_map)

# Export map plot (included .png and .pdf format)
# Uncomment lines below to save
# ggplot2::ggsave("../3_Output/6_other/patients_sample_map.pdf", plot = p_map, height = 6, width = 12)
# ggplot2::ggsave("../3_Output/6_other/patients_sample_map.png", plot = p_map, height = 6, width = 12)

# ====================================================================
# Figure S2: Distribution of study articles and psychiatric patient sample sizes
# ====================================================================

# Statistical Data Preparation
# Reorder factor levels by frequency using forcats
df_stats <- Sample_Info %>% 
  dplyr::filter(!base::is.na(First_Author)) %>%
  dplyr::mutate(Disease_name_G2 = forcats::fct_infreq(Disease_name_G2))

# Descriptive statistics in console
base::print(base::table(df_stats$Disease_name_G2))
total_studies <- base::sum(base::table(df_stats$First_Author))
base::print(base::paste("Total number of studies analyzed:", total_studies))

# Generate bar and jitter plot
p_bar <- ggplot2::ggplot(df_stats, ggplot2::aes(x = Disease_name_G2)) +
  # Bar layer: counts of studies
  ggplot2::geom_bar(
    ggplot2::aes(y = ggplot2::after_stat(count)), 
    fill = "gray70", color = "black", linewidth = 0.3, width = 0.5, alpha = 0.8
  ) +
  # Jitter layer: individual patient sample sizes
  ggplot2::geom_jitter(
    ggplot2::aes(y = Sample_patients, color = Disease_Diagnostic_Criteria_G2),
    position = ggplot2::position_jitterdodge(dodge.width = 0.5, jitter.width = 0.15),
    size = 2.5, alpha = 0.7
  ) +
  # Scales and Dual-Axis (using scales package for breaks)
  ggplot2::scale_y_continuous(
    limits = c(0, 60), expand = c(0, 0),
    breaks = scales::breaks_extended(n = 6),
    sec.axis = ggplot2::dup_axis(name = "Number of studies")
  ) +
  ggplot2::scale_color_brewer(palette = "Paired") +
  # Finalize theme and labels
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 10),
    axis.line.y.right = ggplot2::element_line(color = "black"),
    legend.position = "right",
    legend.background = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.2)
  ) +
  ggplot2::labs(
    x = "Disorder Type", 
    y = "Number of Psychiatric Patients", 
    color = "Diagnostic Criteria"
  )

base::print(p_bar)

# Export map plot (included .png and .pdf format)
# Uncomment lines below to save
# ggplot2::ggsave("../3_Output/6_other/bar_plot.pdf", plot = p_map, height = 6, width = 12)
# ggplot2::ggsave("../3_Output/6_other/bar_plot.png", plot = p_map, height = 6, width = 12)
