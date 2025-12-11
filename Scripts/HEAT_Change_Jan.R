################################################################################
################################################################################
#########################  FWF - Site Prep Study   #############################
#########################    Nov to Jan HEAT       #############################
#########################  University of Florida   #############################
#########################     Gage LaPierre        #############################
#########################     2024 - 2025          #############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014")

#########################      Installs Packages   #############################

list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################      Loads Packages      ############################
library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)

##########################      Read in  Data      #############################
Data = read.csv("Data/Post-Treatment Data.csv")

# Transform data from wide to long format
df_long <- Data %>%
  pivot_longer(
    cols = c('X0m', 'X3m', 'X6m', 'X9m', 'X12m', 'X15m', 'X18m',
             'X21m', 'X24m', 'X27m'),
    names_to = 'Quadrat_ID',
    values_to = 'Coverage'
  ) %>%
  # Remove 'X' and 'm' from Quadrat_ID and create Quadrat_ID_Numeric
  mutate(Quadrat_ID_Numeric = str_replace(Quadrat_ID, 'X', '') %>%
           str_replace('m', '')) %>%
  # Create the unique ID, now including 'Month'
  mutate(Unique_ID = paste(Block, Plot, Quadrat_ID_Numeric, sep = '-'))

# Filter out *only* blank entries (empty strings) in Coverage
# We will handle "NA" and other NA values by converting them to 0
df_cleaned <- df_long %>%
  filter(Coverage != "")

# Convert Coverage to numeric. This will turn any non-numeric strings (including "NA" if still present) into actual R `NA` values.
df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)

# Replace all NA values in Coverage with 0
df_cleaned <- df_cleaned %>%
  mutate(Coverage = replace_na(Coverage, 0))

# Reclassify coverage data (CV) from 1-10 scale to percent scale
# Explicitly ensure the output of case_when is numeric
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

# Orders years and treatments so that they display in same sequence in graphs #
df_cleaned$Month = factor(df_cleaned$Month, levels=c('November','Jan'))

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand',
                                       'Preemergent','Soil Inversion'))

#################### Species abundances ########################################
# Creates and joins  data month to make long data format #
Two_Abundance <- df_cleaned[which(df_cleaned$Month == "November"),]
Three_Abundance <- df_cleaned[which(df_cleaned$Month == "Jan"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance,
                         by = c('Unique_ID', "Treatment", 
                                'Species', 'Block', 'Plot')) # Added Block and Plot to join by
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0,
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0,
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>%
  dplyr::select(Block, Plot, Unique_ID, Treatment, Species,
                Coverage.x, Coverage.y) %>%
  group_by(Block, Plot, Unique_ID, Treatment, Species) %>%
  mutate(Change_abundance = Coverage.y - Coverage.x)


# AGGREGATE TO PLOT LEVEL
Plot_Change_Abundance <- Change_Abundance %>%
  group_by(Block, Plot, Treatment, Species) %>%
  summarise(Mean_Change_abundance = mean(Change_abundance, na.rm = TRUE), 
            .groups = 'drop') %>%
  # Create a unique ID for plots
  mutate(Plot_ID = paste(Block, Plot, sep = '-'))

Plot_Change_Abundance_H = filter(Plot_Change_Abundance, Mean_Change_abundance >= 0)
Plot_Change_Abundance_L = filter(Plot_Change_Abundance, Mean_Change_abundance <= 0)

Plot_Change_Abundance = full_join(Plot_Change_Abundance_H, Plot_Change_Abundance_L)

# Filter out species with less than 5% absolute cover change
Plot_Change_Abundance <- Plot_Change_Abundance %>%
  filter(abs(Mean_Change_abundance) >= 5)

Treat_Plot = ungroup(Plot_Change_Abundance) %>%
  dplyr::select(Plot_ID, Treatment) %>%
  group_by(Treatment, Plot_ID) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Plot_ID')

Data_change_plot <- ungroup(Plot_Change_Abundance) %>%
  dplyr::select(Plot_ID, Species, Mean_Change_abundance) %>%
  as.data.frame() %>%
  matrify()

Data_change_plot_matrix <- as.matrix(Data_change_plot)
rownames(Data_change_plot_matrix) <- Treat_Plot %>% rownames()


# Define the colors for each treatment (for the annotation bar)
annotation_colors = list(
  Treatment = c(Control = "yellow", `Fill Sand` = "orange", 
                Preemergent = "#66CC00", `Soil Inversion` = "#CC66CC")
)

# Get the minimum and maximum values from your matrix
min_val <- min(Data_change_plot_matrix)
max_val <- max(Data_change_plot_matrix)

# Create a sequence of breaks for the color scale, centered at 0
breaks <- c(seq(min_val, 0, length.out = 50),
            seq(0.01, max_val, length.out = 50))

# Define two color palettes: one for negative values and one for positive values
neg_colors <- colorRampPalette(c("blue", "lightblue", "white"))(length(breaks[breaks <= 0]))
pos_colors <- colorRampPalette(c("White", "pink", "red"))(length(breaks[breaks > 0]))

# Combine the two palettes into a single vector of colors for the heatmap
my_heatmap_colors <- c(neg_colors, pos_colors)

# Update your pheatmap function call with the new color and breaks arguments
speciesHEAT_plot = pheatmap(Data_change_plot_matrix, show_rownames = T,
                            cluster_cols = F, cluster_rows = F,
                            annotation_row = Treat_Plot, fontsize = 10,
                            border_color = "black", display_numbers = FALSE,
                            cellheight = 10, cellwidth = 28,
                            # Use the new variable for annotation colors
                            annotation_colors = annotation_colors,
                            # Use the new variable for heatmap cell colors
                            color = my_heatmap_colors,
                            breaks = breaks)

speciesHEAT_plot

png(file = "Figures/Heat_Change_Jan.png",
    units="cm", width=30, height=15, res=100)
speciesHEAT_plot
dev.off()
