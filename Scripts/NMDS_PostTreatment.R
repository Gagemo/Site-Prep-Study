################################################################################
################################################################################
#########################   Site Prep Study       ##############################
#########################    NMDS - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2024 - 2025         ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggrepel","ggsignif", "multcompView", "ggpubr", 
                      "rstatix", 'rmarkdown', "labdsv", "pairwiseAdonis", 
                      "devtools", "knitr", "tables", "openxlsx", "labdsv")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(extrafont)
#font_import()
loadfonts(device = "win")
library(tidyverse)
library(vegan)
library(agricolae)
library(rmarkdown)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(rstatix)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(vegan)
library(labdsv)
library(devtools)
library(knitr)
library(tables)
library(openxlsx)
library(labdsv)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Load necessary libraries
library(tidyverse)
library(labdsv) # Ensure labdsv is loaded for matrify

########################## Read in Data ########################################
Data = read.csv("Data/Post-Treatment Data.csv")

# Filter by Month
Data_April = Data %>% filter(Month == "April")
Data_Jan = Data %>% filter(Month == "Jan")
Data_Nov = Data %>% filter(Month == "November")

################################################################################
################################################################################
################################### April ######################################
################################################################################
################################################################################

# Transform data from wide to long format
df_long <- Data_April %>%
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
  mutate(Unique_ID = paste(Month, Block, Plot, Quadrat_ID_Numeric, sep = '-'))

# Filter out blank and "NA" entries in Coverage
df_cleaned <- df_long %>%
  filter(Coverage != "", Coverage != "NA")

# Convert Coverage to numeric.
df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)

# Remove any NAs introduced by the above numeric coercion.
df_cleaned <- df_cleaned %>% drop_na(Coverage)

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

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

# Create species pivot table
df_for_matrify <- df_cleaned %>%
  dplyr::select(Unique_ID, Species, Coverage) %>%
  mutate(Coverage = as.numeric(Coverage))

# Convert the tibble to a base R data.frame before calling matrify()
df_for_matrify_base <- as.data.frame(df_for_matrify)

# Now, matrify() should work correctly with the standard data.frame
Spp = matrify(df_for_matrify_base)

Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(df_cleaned, Unique_ID, Month, Block, Plot, Treatment) %>% 
  dplyr::summarize()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Veg_Spp, distance = 'bray', k=2)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,
                  Plot = Treat$Plot)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)
envfit_result

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/NMDS_FitValues_April.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp_filtered <- Spp[, -absent_species]

# Perform Indicator Species Analysis
indicator_results <- indval(Spp_filtered, Treat$Treatment)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_BH = indval_df$BH,
  Indicator_Value_BM = indval_df$BM,
  Indicator_Value_LH = indval_df$LH,
  Indicator_Value_LM = indval_df$LM,
  Indicator_Value_W = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, 
      significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/indicator_species_analysis_results_April.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################

simper_result <- simper(Spp, Treat$Treatment)
simper_summary <- summary(simper_result)

# Extract the first pairwise comparison (adjust as needed)
simper_result <- simper(Spp, Treat$Treatment, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/SIMPER_April.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

# NMDS Graphs
NMDS_graph_April = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
 # geom_text(data = species.scores, aes(x = MDS1, y = MDS2, label = species))+
  
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 
                              'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  annotate("text", x = -.9, y = .9,
           label = paste0("Stress: ", format(MDS$stress, digits = 2)),
           hjust = 0, size = 8) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(
    #   fill = guide_legend(nrow = 1), # Removed fill
    shape = guide_legend(nrow = 1)
  ) +
  labs(x = "MDS1", y = "MDS2", title = "6-month",
       color = "", # Added color to the labs
       fill = "",
       shape = "") # added fill to the labs
NMDS_graph_April

NMDS_graph_April_Spp = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2, label = species),
                  fontface = "bold",
                  force = 2, size = 6,
                  max.overlaps = Inf) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  labs(x = "MDS1", y = "MDS2", title = "6-month",) +
  scale_y_continuous(expand = expansion(mult = 0.1))
NMDS_graph_April_Spp

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
write.csv.tabular(adon.results, "Figures/adonis_April.csv")
pairwise.adonis <-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

#save tables
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, 
            startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], 
              startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/pairwise_adonis_same_sheet_April.xlsx", 
             overwrite = TRUE)

################################################################################
################################################################################
################################### Jan ########################################
################################################################################
################################################################################

# Transform data from wide to long format
df_long <- Data_Jan %>%
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
  mutate(Unique_ID = paste(Month, Block, Plot, Quadrat_ID_Numeric, sep = '-'))

# Filter out blank and "NA" entries in Coverage
df_cleaned <- df_long %>%
  filter(Coverage != "", Coverage != "NA")

# Convert Coverage to numeric.
df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)

# Remove any NAs introduced by the above numeric coercion.
df_cleaned <- df_cleaned %>% drop_na(Coverage)

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

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

# Create species pivot table
df_for_matrify <- df_cleaned %>%
  dplyr::select(Unique_ID, Species, Coverage) %>%
  mutate(Coverage = as.numeric(Coverage))

# Convert the tibble to a base R data.frame before calling matrify()
df_for_matrify_base <- as.data.frame(df_for_matrify)

# Now, matrify() should work correctly with the standard data.frame
Spp = matrify(df_for_matrify_base)

Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(df_cleaned, Unique_ID, Month, Block, Plot, Treatment) %>% 
  dplyr::summarize()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Veg_Spp, distance = 'bray', k=2)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,
                  Plot = Treat$Plot)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)
envfit_result

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/NMDS_FitValues_Jan.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp_filtered <- Spp[, -absent_species]

# Perform Indicator Species Analysis
indicator_results <- indval(Spp_filtered, Treat$Treatment)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_BH = indval_df$BH,
  Indicator_Value_BM = indval_df$BM,
  Indicator_Value_LH = indval_df$LH,
  Indicator_Value_LM = indval_df$LM,
  Indicator_Value_W = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, 
      significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/indicator_species_analysis_results_Jan.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################

simper_result <- simper(Spp, Treat$Treatment)
simper_summary <- summary(simper_result)

# Extract the first pairwise comparison (adjust as needed)
simper_result <- simper(Spp, Treat$Treatment, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/SIMPER_Jan.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

# NMDS Graphs
NMDS_graph_Jan = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  # geom_text(data = species.scores, aes(x = MDS1, y = MDS2, label = species))+
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 
                              'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  annotate("text", x = -1, y = 1,
           label = paste0("Stress: ", format(MDS$stress, digits = 2)),
           hjust = 0, size = 8) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(
    #   fill = guide_legend(nrow = 1), # Removed fill
    shape = guide_legend(nrow = 1)
  ) +
  labs(x = "", y = "MDS2", title = "3-month",
       color = "", # Added color to the labs
       fill = "",
       shape = "") # added fill to the labs
NMDS_graph_Jan

NMDS_graph_Jan_Spp = 
ggplot() +
  geom_text_repel(data = species.scores, 
                  aes(x = MDS1, y = MDS2, label = species),
                  fontface = "bold", force = 2, size = 6,
                  max.overlaps = Inf) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
        ,
        legend.position="bottom") +
  labs(x = "", y = "", title = "3-month",) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1))
NMDS_graph_Jan_Spp

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
write.csv.tabular(adon.results, "Figures/adonis_Jan.csv")
pairwise.adonis <-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

#save tables
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, 
            startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], 
              startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/pairwise_adonis_same_sheet_Jan.xlsx", 
             overwrite = TRUE)

################################################################################
################################################################################
########################## November ############################################
################################################################################
################################################################################

# Transform data from wide to long format
df_long <- Data_Nov %>%
  pivot_longer(
    cols = c('X0m', 'X3m', 'X6m', 'X9m', 'X12m', 'X15m', 'X18m',
             'X21m', 'X24m', 'X27m'),
    names_to = 'Quadrat_ID',
    values_to = 'Coverage'
  ) %>%
  # Remove 'X' and 'm' from Quadrat_ID and create Quadrat_ID_Numeric
  mutate(Quadrat_ID_Numeric = str_replace(Quadrat_ID, 'X', '') %>% str_replace('m', '')) %>%
  # Create the unique ID, now including 'Month'
  mutate(Unique_ID = paste(Month, Block, Plot, Quadrat_ID_Numeric, sep = '-'))

# Filter out blank and "NA" entries in Coverage
df_cleaned <- df_long %>%
  filter(Coverage != "", Coverage != "NA")

# Convert Coverage to numeric. This will introduce NAs for truly non-numeric values.
df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)

# Remove any NAs introduced by the above numeric coercion.
df_cleaned <- df_cleaned %>% drop_na(Coverage)

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
    TRUE ~ Coverage # If already numeric and not one of the reclassified values, keep as is
  )))

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

# Create species pivot table
df_for_matrify <- df_cleaned %>%
  dplyr::select(Unique_ID, Species, Coverage) %>%
  mutate(Coverage = as.numeric(Coverage)) # Redundant but safe final numeric coercion

# Convert the tibble to a base R data.frame before calling matrify()
df_for_matrify_base <- as.data.frame(df_for_matrify)

# Now, matrify() should work correctly with the standard data.frame
Spp = matrify(df_for_matrify_base)

Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(df_cleaned, Unique_ID, Month, Block, Plot, Treatment) %>% 
  dplyr::summarize()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Veg_Spp, distance = 'bray', k=2)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,
                  Plot = Treat$Plot)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)
envfit_result

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/NMDS_FitValues_Nov.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp_filtered <- Spp[, -absent_species]

# Perform Indicator Species Analysis
indicator_results <- indval(Spp_filtered, Treat$Treatment)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_BH = indval_df$BH,
  Indicator_Value_BM = indval_df$BM,
  Indicator_Value_LH = indval_df$LH,
  Indicator_Value_LM = indval_df$LM,
  Indicator_Value_W = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, 
      significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/indicator_species_analysis_results_Nov.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################

simper_result <- simper(Spp, Treat$Treatment)
simper_summary <- summary(simper_result)

# Extract the first pairwise comparison (adjust as needed)
simper_result <- simper(Spp, Treat$Treatment, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/SIMPER_Nov.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

# NMDS Graphs
Nov_NMDS_graph = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  #geom_text(data = NMDS, 
  #          aes(x = MDS.MDS1, y = MDS.MDS2, label = Plot), size = 8) +
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 
                              'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  annotate("text", x = -0.19, y = 0.19,
           label = paste0("Stress: ", format(MDS$stress, digits = 2)),
           hjust = 0, size = 8) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(
    #   fill = guide_legend(nrow = 1), # Removed fill
    shape = guide_legend(nrow = 1)
  ) +
  labs(x = "", y = "MDS2", title = "1-month",
       color = "", # Added color to the labs
       fill = "",
       shape = "") # added fill to the labs
Nov_NMDS_graph

Nov_NMDS_graph_Spp =
ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2, label = species),
                  fontface = "bold",
                  force = 2, size = 6,
                  max.overlaps = Inf) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  labs(x = "", y = "", title = "1-month") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1))
Nov_NMDS_graph_Spp

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
write.csv.tabular(adon.results, "Figures/adonis_Nov.csv")
pairwise.adonis <-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

#save tables
# Create a new workbook

wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, 
            startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], 
              startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/pairwise_adonis_same_sheet_Nov.xlsx", overwrite = TRUE)


################################################################################
################################################################################
NMDS_Combined = ggarrange(Nov_NMDS_graph, Nov_NMDS_graph_Spp,
                          NMDS_graph_Jan , NMDS_graph_Jan_Spp,
                          NMDS_graph_April, NMDS_graph_April_Spp,
                          nrow = 3, ncol = 2, common.legend = TRUE,
                          legend = "bottom")
NMDS_Combined
ggsave("Figures/NMDS_Combined.png", dpi = 1000, width = 18, height = 24)
