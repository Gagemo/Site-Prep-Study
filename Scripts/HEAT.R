################################################################################
################################################################################
#########################  FWF - Site Prep Study   #############################
#########################   HEATMAP + BARE GROUND  #############################
#########################  University of Florida   #############################
#########################      Gage LaPierre       #############################
#########################       2024 - 2025         #############################
################################################################################
################################################################################

######################### Clears Environment & History #########################
rm(list = ls(all = TRUE))
cat("\014")

######################### Installs Packages #####################################
list.of.packages <- c(
  "tidyverse", "vegan", "labdsv",
  "pheatmap", "cowplot", "grid"
)

new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages)) {
  install.packages(new.packages)
}

######################### Loads Packages ########################################
library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)
library(cowplot)
library(grid)

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

################################################################################
# DATA PROCESSING
################################################################################

process_data_for_month <- function(month) {
  
  Data <- read.csv(
    "Data/Post-Treatment Data.csv",
    stringsAsFactors = FALSE,
    check.names = TRUE
  ) %>%
    dplyr::filter(Month == month)
  
  df_long <- Data %>%
    tidyr::pivot_longer(
      cols = c(
        "X0m", "X3m", "X6m", "X9m", "X12m",
        "X15m", "X18m", "X21m", "X24m", "X27m"
      ),
      names_to = "Quadrat_ID",
      values_to = "Coverage"
    ) %>%
    dplyr::mutate(
      Quadrat_ID_Numeric = stringr::str_replace(Quadrat_ID, "X", "") %>%
        stringr::str_replace("m", "")
    )
  
  # Blank entries, ".", NA, and N/A are treated as absence and converted to zero.
  df_cleaned <- df_long %>%
    dplyr::mutate(
      Coverage = trimws(as.character(Coverage)),
      Coverage = dplyr::case_when(
        is.na(Coverage) ~ "0",
        Coverage == "" ~ "0",
        toupper(Coverage) %in% c("NA", "N/A", ".") ~ "0",
        TRUE ~ Coverage
      ),
      Coverage = suppressWarnings(as.numeric(Coverage))
    )
  
  if (any(is.na(df_cleaned$Coverage))) {
    stop(
      "Some Coverage values could not be converted to numeric for month: ",
      month
    )
  }
  
  # Convert cover classes to representative percentage midpoints.
  df_cleaned <- df_cleaned %>%
    dplyr::mutate(
      Coverage = dplyr::case_when(
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
        TRUE ~ NA_real_
      ),
      Treatment = factor(
        Treatment,
        levels = c(
          "Control",
          "Fill Sand",
          "Preemergent",
          "Soil Inversion"
        )
      ),
      Species = trimws(as.character(Species))
    )
  
  if (any(is.na(df_cleaned$Coverage))) {
    stop(
      "Coverage values outside the expected 0-10 scale were found for month: ",
      month
    )
  }
  
  # Average the ten quadrats within each plot.
  df_by_plot <- df_cleaned %>%
    dplyr::filter(Species != "") %>%
    dplyr::group_by(
      Block,
      Plot,
      Treatment,
      Species
    ) %>%
    dplyr::summarise(
      Mean_Coverage = mean(
        Coverage,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Month = month,
      Plot_ID = paste(
        Block,
        Plot,
        sep = "-"
      )
    ) %>%
    dplyr::arrange(
      Treatment,
      Block,
      Plot
    )
  
  return(df_by_plot)
}

######################### Process Sampling Periods ##############################
df_by_plot_Jan <- process_data_for_month("Jan")
df_by_plot_April <- process_data_for_month("April")

all_plot_data <- dplyr::bind_rows(
  df_by_plot_Jan,
  df_by_plot_April
)

################################################################################
# SEPARATE BARE GROUND FROM SPECIES DATA
################################################################################

# This pattern recognizes common versions of the bare-ground label:
# "Bare Ground", "Bareground", "Bare-ground", and differences in capitalization.
bare_ground_pattern <- "^bare[[:space:]_-]*ground$"

bare_ground_data <- all_plot_data %>%
  dplyr::filter(
    stringr::str_detect(
      stringr::str_to_lower(Species),
      bare_ground_pattern
    )
  )

species_data <- all_plot_data %>%
  dplyr::filter(
    !stringr::str_detect(
      stringr::str_to_lower(Species),
      bare_ground_pattern
    )
  )

if (nrow(bare_ground_data) == 0) {
  stop(
    "No bare-ground rows were found. Check the spelling used in the Species column."
  )
}

################################################################################
# RETAIN SPECIES FOR THE HEATMAP
################################################################################

# Retain a species if its mean cover reaches at least 5% in at least one plot
# during either sampling period. The same species set is then used in both
# heatmaps so the panels are directly comparable.

species_to_keep <- species_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Maximum_Plot_Mean = max(
      Mean_Coverage,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(Maximum_Plot_Mean >= 5) %>%
  dplyr::pull(Species)

species_data <- species_data %>%
  dplyr::filter(
    Species %in% species_to_keep
  )

species_data_Jan <- species_data %>%
  dplyr::filter(Month == "Jan")

species_data_April <- species_data %>%
  dplyr::filter(Month == "April")

# Shared color limits are calculated from species cover only.
# Bare ground no longer compresses the species color scale.
min_val <- 0
max_val <- max(
  species_data$Mean_Coverage,
  na.rm = TRUE
)

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

create_heatmap_plot <- function(
    data,
    title,
    species_order,
    show_annotation = TRUE,
    show_plot_labels = TRUE
) {
  
  # Complete the plot-by-species matrix so taxa absent from a plot equal zero.
  plot_metadata <- data %>%
    dplyr::select(
      Plot_ID,
      Block,
      Plot,
      Treatment
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(
      Treatment,
      Block,
      Plot
    )
  
  heatmap_matrix <- data %>%
    dplyr::select(
      Plot_ID,
      Species,
      Mean_Coverage
    ) %>%
    tidyr::complete(
      Plot_ID = plot_metadata$Plot_ID,
      Species = species_order,
      fill = list(
        Mean_Coverage = 0
      )
    ) %>%
    dplyr::mutate(
      Species = factor(
        Species,
        levels = species_order
      )
    ) %>%
    tidyr::pivot_wider(
      names_from = Species,
      values_from = Mean_Coverage,
      values_fill = 0
    ) %>%
    dplyr::right_join(
      plot_metadata %>%
        dplyr::select(Plot_ID),
      by = "Plot_ID"
    )
  
  heatmap_matrix_data <- heatmap_matrix %>%
    dplyr::arrange(
      match(
        Plot_ID,
        plot_metadata$Plot_ID
      )
    ) %>%
    tibble::column_to_rownames("Plot_ID") %>%
    as.matrix()
  
  annotation_data <- plot_metadata %>%
    dplyr::select(
      Plot_ID,
      Treatment
    ) %>%
    tibble::column_to_rownames("Plot_ID")
  
  final_annotation_row <- if (show_annotation) {
    annotation_data
  } else {
    NULL
  }
  
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
    color = colorRampPalette(
      c("white", "red")
    )(50),
    breaks = seq(
      min_val,
      max_val,
      length.out = 51
    ),
    legend = FALSE,
    annotation_legend = FALSE
  )
  
  return(p_plot)
}

# Keep the species columns in the same order in both panels.
# Species are ordered from highest to lowest mean cover across both periods.
species_order <- species_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Overall_Mean = mean(
      Mean_Coverage,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Overall_Mean)
  ) %>%
  dplyr::pull(Species)

speciesHEAT_plot_Jan <- create_heatmap_plot(
  species_data_Jan,
  "3-month",
  species_order = species_order,
  show_annotation = FALSE,
  show_plot_labels = TRUE
)

speciesHEAT_plot_April <- create_heatmap_plot(
  species_data_April,
  "6-month",
  species_order = species_order,
  show_annotation = TRUE,
  show_plot_labels = FALSE
)

######################### Shared Color Legend ###################################

legend_data_plot <- ggplot(
  data.frame(
    Mean_Coverage = seq(
      min_val,
      max_val,
      length.out = 100
    ),
    x = 1
  ),
  aes(
    x = x,
    y = Mean_Coverage,
    fill = Mean_Coverage
  )
) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorRampPalette(
      c("white", "red")
    )(50),
    limits = c(
      min_val,
      max_val
    ),
    name = "Mean species\ncover (%)"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )

color_legend_grob <- cowplot::get_legend(
  legend_data_plot
)

######################### Treatment Legend ######################################

treatment_data <- data.frame(
  Treatment = factor(
    c(
      "Control",
      "Fill Sand",
      "Preemergent",
      "Soil Inversion"
    ),
    levels = c(
      "Control",
      "Fill Sand",
      "Preemergent",
      "Soil Inversion"
    )
  ),
  Dummy = 1
)

treatment_legend_plot <- ggplot(
  treatment_data,
  aes(
    x = Dummy,
    y = Treatment,
    fill = Treatment
  )
) +
  geom_tile() +
  scale_fill_manual(
    values = annotation_colors$Treatment,
    name = "Treatment"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )

treatment_legend_grob <- cowplot::get_legend(
  treatment_legend_plot
)

legends_plot <- cowplot::plot_grid(
  treatment_legend_grob,
  color_legend_grob,
  ncol = 1,
  align = "v",
  rel_heights = c(
    0.6,
    0.5
  )
)

######################### Combine Species Heatmaps ##############################

heat_Jan <- cowplot::ggdraw() +
  cowplot::draw_grob(
    speciesHEAT_plot_Jan$gtable
  )

heat_April <- cowplot::ggdraw() +
  cowplot::draw_grob(
    speciesHEAT_plot_April$gtable
  )

combined_species_heatmap <- cowplot::plot_grid(
  heat_Jan,
  heat_April,
  legends_plot,
  nrow = 1,
  rel_widths = c(
    0.9,
    0.9,
    0.6
  ),
  align = "h",
  axis = "l"
)

ggsave(
  "Figures/Combined_Species_Heatmaps.tiff",
  plot = combined_species_heatmap,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

################################################################################
# SEPARATE BARE-GROUND FIGURE
################################################################################

bare_ground_data <- bare_ground_data %>%
  dplyr::mutate(
    Sampling_Period = factor(
      Month,
      levels = c(
        "Jan",
        "April"
      ),
      labels = c(
        "3-month",
        "6-month"
      )
    )
  )

bare_ground_plot <- ggplot(
  bare_ground_data,
  aes(
    x = Treatment,
    y = Mean_Coverage,
    fill = Treatment
  )
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  geom_jitter(
    color = "black",
    size = 2.5,
    alpha = 0.75,
    width = 0.15,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Sampling_Period,
    nrow = 1
  ) +
  scale_fill_manual(
    values = annotation_colors$Treatment
  ) +
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    x = "",
    y = "Mean bare-ground cover (%)"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      size = 16
    ),
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 13,
      face = "bold",
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 13,
      face = "bold"
    )
  )

print(bare_ground_plot)

ggsave(
  "Figures/Bare_Ground_Cover.tiff",
  plot = bare_ground_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)