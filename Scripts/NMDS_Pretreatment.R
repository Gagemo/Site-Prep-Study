################################################################################
################################################################################
#########################    Site Prep Study      ##############################
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
                      "devtools", "knitr", "tables", "openxlsx", "labdsv", 
                      'ggrepel')
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
library(ggrepel)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Load necessary libraries
library(tidyverse)
library(labdsv) # Ensure labdsv is loaded for matrify

########################## Read in Data ########################################
Data = read.csv("Data/Pretreatment Data.csv")

################################################################################
################################################################################
################################### June  ######################################
################################################################################
################################################################################

# Filter out blank and "NA" entries in Coverage
df_cleaned <- Data %>%
  filter(Coverage != "", Coverage != "NA")

# Create the unique ID, now including 'Month'
df_cleaned <- df_cleaned %>% mutate(Unique_ID = paste(Block, 
                                       Plot, Quadrat, sep = '-'))

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
                              levels=c('Control', 'Fill', 
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
MDS = metaMDS(Veg_Spp, distance = 'bray', k=3)
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
                  Block = Treat$Block, Plot = Treat$Plot, Unique_ID = Treat$Unique_ID)

################################################################################
#############################NMDS Graphs########################################
################################################################################

# NMDS Graphs
NMDS_graph = ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
 # geom_text_repel(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, label = Unique_ID),
  #                size = 3) + # Replace geom_text with geom_text_repel
 # geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2, label = species),
 #                 size = 3) + # Replace geom_text with geom_text_repel
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(shape = guide_legend(nrow = 1)) +
  labs(x = "MDS1", y = "MDS2", color = "Treatment", fill = "Treatment")
NMDS_graph

ggsave("Figures/pretreat_NMDS_graph.tiff", dpi = 100)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
write.csv.tabular(adon.results, "Tables/Pretreat_adonis.csv")
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
saveWorkbook(wb, "Tables/Pretreat_pairwise_adonis_same_sheet.xlsx", 
             overwrite = TRUE)