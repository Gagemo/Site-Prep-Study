################################################################################
################################################################################
#########################  FWF - Site Prep Study   #############################
#########################            HEAT          #############################
#########################  University of Florida   #############################
#########################      Gage LaPierre       #############################
#########################       2024 - 2025        #############################
################################################################################
################################################################################

# Clears Environment & History
rm(list = ls(all = TRUE))
cat("\014")

# Install and load necessary packages
list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap", "cowplot", "grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)
library(cowplot)
library(grid)

################################################################################
# DATA PROCESSING
################################################################################

process_data_for_month <- function(month) {
  
  Data <- read.csv("Data/Post-Treatment Data.csv", stringsAsFactors = FALSE) %>%
    dplyr::filter(Month == month)
  
  df_long <- Data %>%
    tidyr::pivot_longer(
      cols = c("X0m", "X3m", "X6m", "X9m", "X12m", "X15m", "X18m",
               "X21m", "X24m", "X27m"),
      names_to = "Quadrat_ID",
      values_to = "Coverage"
    ) %>%
    dplyr::mutate(
      Quadrat_ID_Numeric = stringr::str_replace(Quadrat_ID, "X", "") %>%
        stringr::str_replace("m", "")
    )
  
  # Robust cleaning/parsing for Coverage
  # Treat "." as missing (common data-entry placeholder), then parse and replace NA with 0
  df_cleaned <- df_long %>%
    dplyr::mutate(Coverage = trimws(Coverage)) %>%
    dplyr::mutate(Coverage = dplyr::na_if(Coverage, "")) %>%
    dplyr::mutate(Coverage = dplyr::na_if(Coverage, ".")) %>%
    dplyr::mutate(Coverage = readr::parse_number(Coverage)) %>%
    dplyr::mutate(Coverage = tidyr::replace_na(Coverage, 0))
  
  # Convert Braun-Blanquet-ish classes to midpoints
  df_cleaned <- df_cleaned %>%
    dplyr::mutate(
      Coverage = as.numeric(dplyr::case_when(
        Coverage == 10 ~ 97.5,
        Coverage == 1  ~ 0.1,
        Coverage == 2  ~ 0.5,
        Coverage == 3  ~ 1.5,
        Coverage == 4  ~ 3.5,
        Coverage == 5  ~ 7.5,
        Coverage == 6  ~ 17.5,
        Coverage == 7  ~ 37.5,
        Coverage == 8  ~ 62.5,
        Coverage == 9  ~ 85,
        Coverage == 0  ~ 0,
        TRUE ~ Coverage
      ))
    )
  
  df_cleaned$Treatment <- factor(
    df_cleaned$Treatment,
    levels = c("Control", "Fill Sand", "Preemergent", "Soil Inversion")
  )
  
  df_by_plot <- df_cleaned %>%
    dplyr::filter(Species != "") %>%
    dplyr::group_by(Block, Plot, Treatment, Species) %>%
    dplyr::summarise(Mean_Coverage = mean(Coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Treatment)
  
  # Keep only species with >= 5% mean cover
  df_by_plot <- df_by_plot %>%
    dplyr::filter(abs(Mean_Coverage) >= 5)
  
  return(df_by_plot)
}

# Process data for both months
df_by_plot_Jan <- process_data_for_month("Jan")
df_by_plot_April <- process_data_for_month("April")

# Combine data to find a universal max/min for the color scale
all_data <- dplyr::bind_rows(df_by_plot_Jan, df_by_plot_April)
min_val <- min(all_data$Mean_Coverage, na.rm = TRUE)
max_val <- max(all_data$Mean_Coverage, na.rm = TRUE)

################################################################################
# HEATMAP + LEGENDS
################################################################################

annotation_colors <- list(
  Treatment = c(
    Control = "yellow",
    `Fill Sand` = "orange",
    Preemergent = "#66CC00",
    `Soil Inversion` = "#CC66CC"
  )
)

create_heatmap_plot <- function(data, title, show_annotation = TRUE, show_plot_labels = TRUE) {
  
  heatmap_matrix <- data %>%
    dplyr::mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
    dplyr::select(Plot_ID, Species, Mean_Coverage) %>%
    tidyr::pivot_wider(
      names_from = Species,
      values_from = Mean_Coverage,
      values_fill = 0
    )
  
  heatmap_matrix_data <- heatmap_matrix %>%
    tibble::column_to_rownames("Plot_ID") %>%
    as.matrix()
  
  annotation_data <- data %>%
    dplyr::mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
    dplyr::select(Plot_ID, Treatment) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("Plot_ID")
  
  final_annotation_row <- if (show_annotation) annotation_data else NULL
  
  p_plot <- pheatmap::pheatmap(
    heatmap_matrix_data,
    show_rownames = show_plot_labels,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = final_annotation_row,
    annotation_colors = annotation_colors,
    fontsize = 12,
    border_color = "black",
    display_numbers = FALSE,
    cellheight = 15,
    cellwidth = 28,
    main = title,
    color = colorRampPalette(c("white", "red"))(50),
    legend = FALSE,
    annotation_legend = FALSE
  )
  
  return(p_plot)
}

# Heatmaps
speciesHEAT_plot_Jan <- create_heatmap_plot(df_by_plot_Jan, "3-month",
                                            show_annotation = FALSE, show_plot_labels = TRUE)

speciesHEAT_plot_April <- create_heatmap_plot(df_by_plot_April, "6-month",
                                              show_annotation = TRUE, show_plot_labels = FALSE)

# Shared color legend (dummy ggplot)
legend_data_plot <- ggplot(all_data, aes(x = 1, y = Mean_Coverage, fill = Mean_Coverage)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorRampPalette(c("white", "red"))(50),
    limits = c(min_val, max_val),
    name = "Mean coverage (%)"
  ) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))

color_legend_grob <- cowplot::get_legend(legend_data_plot)

# Treatment legend (dummy ggplot)
treatment_data <- data.frame(
  Treatment = c("Control", "Fill Sand", "Preemergent", "Soil Inversion"),
  Dummy = 1
)

treatment_legend_plot <- ggplot(treatment_data, aes(x = Dummy, y = Treatment, fill = Treatment)) +
  geom_tile() +
  scale_fill_manual(values = annotation_colors$Treatment, name = "Treatment") +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))

treatment_legend_grob <- cowplot::get_legend(treatment_legend_plot)

# IMPORTANT: build legends_plot (this is what your run was missing)
legends_plot <- cowplot::plot_grid(
  treatment_legend_grob,
  color_legend_grob,
  ncol = 1,
  align = "v",
  rel_heights = c(0.6, 0.5)
)

# Convert pheatmap gtables into cowplot-friendly panels
heat_Jan <- cowplot::ggdraw() + cowplot::draw_grob(speciesHEAT_plot_Jan$gtable)
heat_April <- cowplot::ggdraw() + cowplot::draw_grob(speciesHEAT_plot_April$gtable)

# Combine
combined_plot <- cowplot::plot_grid(
  heat_Jan,
  heat_April,
  legends_plot,
  nrow = 1,
  rel_widths = c(0.9, 0.9, 0.6),
  align = "h",
  axis = "l"
)

# Save the final combined plot
ggsave(
  "Figures/Combined_Heatmaps.tiff",
  combined_plot,
  width = 14,
  height = 7,
  units = "in",
  dpi = 100
)
