################################################################################
################################################################################
#########################  FWF - Site Prep Study   #############################
#########################     November HEAT        #############################
#########################  University of Florida   #############################
#########################     Gage LaPierre        #############################
#########################     2024 - 2025          #############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014")

# Install and load packages (same as before)
list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)

# Read and clean data (same as before)
Data = read.csv("Data/Post-Treatment Data.csv")
Data = filter(Data, Month == "November")

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

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand',
                                       'Preemergent','Plow'))

# --- NEW CODE ADDED HERE ---

# 1. Aggregate data by Plot and Species (same as before)
df_by_plot <- df_cleaned %>%
  group_by(Block, Plot, Treatment, Species) %>%
  summarise(Mean_Coverage = mean(Coverage, na.rm = TRUE), .groups = 'drop') %>%
  # 2. Add the new step: Order the data frame by the Treatment factor
  # This is the key change that sets the row order for the heatmap
  arrange(Treatment)

# 3. Reshape the aggregated and ordered data for the heatmap
heatmap_matrix <- df_by_plot %>%
  mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
  select(Plot_ID, Species, Mean_Coverage) %>%
  pivot_wider(names_from = Species, values_from = Mean_Coverage, values_fill = 0)

# 4. Separate the `Plot_ID` column and make it the row names
heatmap_matrix_data <- heatmap_matrix %>%
  column_to_rownames("Plot_ID") %>%
  as.matrix()

# 5. Create the annotation data frame, making sure the order matches the matrix
annotation_data <- df_by_plot %>%
  mutate(Plot_ID = paste(Block, Plot, sep = "-")) %>%
  select(Plot_ID, Treatment) %>%
  distinct() %>%
  column_to_rownames("Plot_ID")

# 6. Now, run pheatmap with the new matrix and annotation data
speciesHEAT_plot = pheatmap(
  heatmap_matrix_data,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE, # This is crucial. We turn off clustering to maintain the order we set
  annotation_row = annotation_data,
  fontsize = 8,
  border_color = "black",
  display_numbers = FALSE,
  cellheight = 15,
  cellwidth = 28,
  color = colorRampPalette(c("white","blue", "orange", "red"))(50)
)

# Save the plot
png(file = "Figures/Heat_Nov.png",
    units="cm", width=20, height=20, res=300)
speciesHEAT_plot
dev.off()


