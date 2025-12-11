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
rm(list=ls(all=TRUE))
cat("\014")

# Install and load necessary packages
list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap", "cowplot", "grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)
library(cowplot)
library(grid)

# Define a function to process the data for a given month
process_data_for_month <- function(month) {
  Data <- read.csv("Data/Post-Treatment Data.csv") %>%
    filter(Month == month)
  
  df_long <- Data %>%
    pivot_longer(
      cols = c('X0m', 'X3m', 'X6m', 'X9m', 'X12m', 'X15m', 'X18m',
               'X21m', 'X24m', 'X27m'),
      names_to = 'Quadrat_ID',
      values_to = 'Coverage'
    ) %>%
    mutate(Quadrat_ID_Numeric = str_replace(Quadrat_ID, 'X', '') %>%
             str_replace('m', ''))
  
  df_cleaned <- df_long %>%
    filter(Coverage != "")
  
  df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)
  df_cleaned <- df_cleaned %>%
    mutate(Coverage = replace_na(Coverage, 0))
  
  df_cleaned <- df_cleaned %>%
    mutate(Coverage = as.numeric(case_when(
      Coverage == 10 ~ 97.5,
      Coverage == 1 ~ 0.1,
      Coverage == 2 ~ 0.5,
      Coverage == 3 ~ 1.5,
      Coverage == 4 ~ 3.5,
      Coverage == 5 ~ 7.5,
      Coverage == 6 ~ 17.5,
      Coverage == 7 ~ 37.5,
      Coverage == 8 ~ 62.5,
      Coverage == 9 ~ 85,
      Coverage == 0 ~ 0,
      TRUE ~ Coverage)))
  
  df_cleaned$Treatment <- factor(df_cleaned$Treatment,
                                 levels=c('Control', 'Fill Sand', 'Preemergent',
                                          'Soil Inversion'))
  
  df_by_plot <- df_cleaned %>%
    filter(Species != "") %>%
    group_by(Block, Plot, Treatment, Species) %>%
    summarise(Mean_Coverage = mean(Coverage, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Treatment)
  
  df_by_plot <- df_by_plot %>%
    filter(abs(Mean_Coverage) >= 5)
  
  return(df_by_plot)
}

# Process data for both months
df_by_plot_Jan <- process_data_for_month("Jan")
df_by_plot_April <- process_data_for_month("April")

# Combine data to find a universal max/min for the color scale
all_data <- bind_rows(df_by_plot_Jan, df_by_plot_April)
min_val <- min(all_data$Mean_Coverage)
max_val <- max(all_data$Mean_Coverage)

# Define the colors for each treatment (for the annotation bar)
annotation_colors = list(
  Treatment = c(Control = "yellow", `Fill Sand` = "orange",  
                Preemergent = "#66CC00", `Soil Inversion` = "#CC66CC")
)

# Function to create heatmap data and plot
# This function generates the individual pheatmap objects without any legends.
# Added optional arguments for annotation_row and show_rownames
create_heatmap_plot <- function(data, title, show_annotation = TRUE, show_plot_labels = TRUE) {
  # Reshape data
  heatmap_matrix <- data %>%
    mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
    select(Plot_ID, Species, Mean_Coverage) %>%
    pivot_wider(names_from = Species, values_from = Mean_Coverage, values_fill = 0)
  
  # Prepare matrix and annotation
  heatmap_matrix_data <- heatmap_matrix %>%
    column_to_rownames("Plot_ID") %>%
    as.matrix()
  
  annotation_data <- data %>%
    mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
    select(Plot_ID, Treatment) %>%
    distinct() %>%
    column_to_rownames("Plot_ID")
  
  # Conditionally set annotation_row
  final_annotation_row <- if (show_annotation) annotation_data else NA
  
  # Create the pheatmap object
  p_plot <- pheatmap(
    heatmap_matrix_data,
    show_rownames = show_plot_labels, # Use the new parameter
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = final_annotation_row, # Use the conditional annotation
    annotation_colors = annotation_colors,
    fontsize = 12,
    border_color = "black",
    display_numbers = FALSE,
    cellheight = 15,
    cellwidth = 28,
    main = title,
    color = colorRampPalette(c("white","red"))(50),
    legend = FALSE, # Remove the main color legend
    annotation_legend = FALSE # Remove the treatment annotation legend
  )
  
  return(p_plot)
}

# Create individual pheatmap objects
# For Jan (3-month): show annotation, show plot labels
speciesHEAT_plot_Jan <- create_heatmap_plot(df_by_plot_Jan, "3-month", show_annotation = FALSE, show_plot_labels = TRUE)
# For April (6-month): hide annotation, hide plot labels
speciesHEAT_plot_April <- create_heatmap_plot(df_by_plot_April, "6-month", show_annotation = TRUE, show_plot_labels = FALSE)

# Create a dummy ggplot object to extract the shared color legend
legend_data_plot <- ggplot(all_data, aes(x = 1, y = Mean_Coverage, fill = Mean_Coverage)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorRampPalette(c("white", "red"))(50),
    limits = c(min_val, max_val),
    name = "Mean coverage (%)"
  ) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))

# Extract the color legend from the dummy plot
color_legend_grob <- get_legend(legend_data_plot)

# Create a dummy ggplot object to get the treatment legend
treatment_data <- data.frame(
  Treatment = c('Control', 'Fill Sand', 'Preemergent', 'Soil Inversion'),
  Dummy = 1
)

# Create the ggplot plot for the treatment legend
treatment_legend_plot <- ggplot(treatment_data, aes(x = Dummy, y = Treatment, fill = Treatment)) +
  geom_tile() +
  scale_fill_manual(
    values = annotation_colors$Treatment,
    name = "Treatment"
  ) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))

# Extract the treatment legend from the ggplot object
treatment_legend_grob <- get_legend(treatment_legend_plot)

# Convert pheatmap objects to grobs for cowplot
grob_Jan <- grid::grid.grabExpr(grid::grid.draw(speciesHEAT_plot_Jan$gtable))
grob_April <- grid::grid.grabExpr(grid::grid.draw(speciesHEAT_plot_April$gtable))

# Arrange the legends one above the other in a separate plot_grid
# Use rel_heights to control the vertical spacing.
legends_plot <- cowplot::plot_grid(
  treatment_legend_grob,
  color_legend_grob,
  ncol = 1,
  align = "v",
  rel_heights = c(0.6, 0.5)
  )
legends_plot

# Arrange the heatmaps and the stacked legends side-by-side
combined_plot <- cowplot::plot_grid(
  grob_Jan,
  grob_April,
  legends_plot,
  nrow = 1,
  rel_widths = c(0.9, 0.3, 0.6),
  align = "h",
  axis = "l"
)

# Save the final combined plot
ggsave("Figures/Combined_Heatmaps.png", combined_plot, width = 14, height = 7, 
       units = "in", dpi = 600)

