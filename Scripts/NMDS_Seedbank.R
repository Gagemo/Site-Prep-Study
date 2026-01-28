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

#########################      Installs & Loads Packages    ####################
# A more streamlined way to check and install packages
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont",
                      "ggrepel","ggsignif", "multcompView", "ggpubr",
                      "rstatix", "rmarkdown", "labdsv", "pairwiseAdonis",
                      "devtools", "knitr", "tables", "openxlsx")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load all packages
invisible(lapply(list.of.packages, library, character.only = TRUE))

# Install pairwiseAdonis from GitHub if not already installed.
# This might be slow, you only need to do this once.
if (!"pairwiseAdonis" %in% installed.packages()) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)

########################## Read in Data ########################################
# Reading the Seedbank.csv directly
Data = read.csv("Data/Seedbank.csv")

# Create a unique ID for each sample based on Month, Block, Plot, and Sample
Data <- Data %>%
  mutate(Unique_ID = paste(Month, Block, Plot, Sample, sep = '-'))

# Convert Treatment to a factor with specified levels
Data$Treatment = factor(Data$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

# Convert Month to a factor with specified levels
Data$Month = factor(Data$Month,
                        levels=c('Pre-Treatment', 'Post-Treatment'))

Pre_Data = Data %>% filter(Month == "Pre-Treatment")
Post_Data = Data %>% filter(Month == "Post-Treatment")

################################################################################
###################### Pre-Treatment Data ######################################
################################################################################

# Separate the metadata (Treatments) from the species abundance data.
# The species data starts from the column after 'Treatment'
# We will use this to get the row names for our species matrix.
metadata <- Pre_Data %>%
  dplyr::select(Unique_ID, Month, Block, Plot, Sample, Treatment)

# The species abundance data are in the columns from 'Ambrosia artemisiifolia' onwards.
# I'll select all columns from the 6th column onwards.
species_matrix <- Pre_Data %>%
  dplyr::select(Ambrosia.artemisiifolia:Triodanis.perfoliata)

# Check that the dimensions match
if(nrow(metadata) != nrow(species_matrix)) {
  stop("Metadata and species data do not have the same number of rows!")
}

# The NMDS function works best with the metadata separated and the species matrix having rownames
# Let's set the row names for the species_matrix to be the Unique_ID from the metadata.
rownames(species_matrix) <- metadata$Unique_ID

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(metadata, Unique_ID, Block, Plot, Treatment) %>% 
  dplyr::summarize()

######################### NMDS Analysis ########################################
# Calculate Bray-Curtis dissimilarity matrix. This is a common choice for
# community ecology data.
Veg_Spp = vegdist(species_matrix, method = 'bray')

# Perform NMDS ordination.
# We'll use k=2 dimensions for easier visualization. You can uncomment the
# NMDS.scree function from your original script to check for optimal dimensions.
# The 'try' function is used to automatically restart the NMDS if it fails to
# converge, which can be common.
MDS = try(metaMDS(species_matrix, distance = 'bray', k = 3, trymax = 100))
if(inherits(MDS, "try-error")) {
  stop("NMDS failed to converge after multiple tries. Consider using different parameters or checking your data.")
}
MDS$stress

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, species_matrix))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment, 
                  Block = Treat$Block, Plot = Treat$Plot, 
                  Unique_ID = Treat$Unique_ID)


############################# NMDS Graphs ######################################
# Create the NMDS plot using ggplot2.
# We will use stat_ellipse to draw confidence ellipses around the treatment groups.
Pre_NMDS_graph = ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  #geom_text_repel(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, label = Unique_ID),
  #                size = 3) + # Replace geom_text with geom_text_repel
  #geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2, label = species),
  #                size = 3) + # Replace geom_text with geom_text_repel
  
  #stat_ellipse(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, group = Month), 
  #            linetype = "solid") +
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  annotate("text", x = -0.5, y = 1,
           label = paste0("Stress: ", format(MDS$stress, digits = 2)),
           hjust = 0, size = 8) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_text(size=25, face="bold", colour = "black"),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="none") +
  labs(x = "MDS1", y = "MDS2", color = "", fill = "", shape = "",
       title = "Pre-treatment")
Pre_NMDS_graph

########################## Adonis Analysis #####################################
# Perform Adonis analysis (Permutational MANOVA) to test for significant
# differences in community structure between treatments.
# We use the Bray-Curtis distance matrix and specify the data frame containing the 'Treatment' variable.
# The data frame provided to `adonis2` must have rows in the same order as the distance matrix.
adon.results <- adonis2(Veg_Spp ~ Treatment, data = metadata, method = "bray", permutations = 999)
print("Adonis Results:")
print(adon.results)

# Save the Adonis results to a CSV file.
write.csv(adon.results, "Tables/pretreat_Seedbank_adonis.csv")

# Perform pairwise Adonis to see which specific treatments are different from each other.
# This is only necessary if the main adonis test is significant.
pairwise.adon.results <- pairwise.adonis2(Veg_Spp ~ Treatment, data = metadata,
                                          permutations = 999, method = "bray")
print("Pairwise Adonis Results:")
print(pairwise.adon.results)

# Save the pairwise Adonis results to a CSV file.
write.csv(pairwise.adon.results, "Tables/pretreat_Seedbank_pairwise_adonis.csv")

################################################################################
###################### Post-Treatment Data ######################################
################################################################################

# Separate the metadata (Treatments) from the species abundance data.
# The species data starts from the column after 'Treatment'
# We will use this to get the row names for our species matrix.
metadata <- Post_Data %>%
  dplyr::select(Unique_ID, Month, Block, Plot, Sample, Treatment)

# The species abundance data are in the columns from 'Ambrosia artemisiifolia' onwards.
# I'll select all columns from the 6th column onwards.
species_matrix <- Post_Data %>%
  dplyr::select(Ambrosia.artemisiifolia:Triodanis.perfoliata)

# Check that the dimensions match
if(nrow(metadata) != nrow(species_matrix)) {
  stop("Metadata and species data do not have the same number of rows!")
}

# The NMDS function works best with the metadata separated and the species matrix having rownames
# Let's set the row names for the species_matrix to be the Unique_ID from the metadata.
rownames(species_matrix) <- metadata$Unique_ID

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(metadata, Unique_ID, Block, Plot, Treatment) %>% 
  dplyr::summarize()

######################### NMDS Analysis ########################################
# Calculate Bray-Curtis dissimilarity matrix. This is a common choice for
# community ecology data.
Veg_Spp = vegdist(species_matrix, method = 'bray')

# Perform NMDS ordination.
# We'll use k=2 dimensions for easier visualization. You can uncomment the
# NMDS.scree function from your original script to check for optimal dimensions.
# The 'try' function is used to automatically restart the NMDS if it fails to
# converge, which can be common.
MDS = try(metaMDS(species_matrix, distance = 'bray', k = 3, trymax = 100))
if(inherits(MDS, "try-error")) {
  stop("NMDS failed to converge after multiple tries. Consider using different parameters or checking your data.")
}
MDS$stress

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, species_matrix))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment, 
                  Block = Treat$Block, Plot = Treat$Plot, 
                  Unique_ID = Treat$Unique_ID)


############################# NMDS Graphs ######################################
# Create the NMDS plot using ggplot2.
# We will use stat_ellipse to draw confidence ellipses around the treatment groups.
Post_NMDS_graph = ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  #geom_text_repel(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, label = Unique_ID),
  #                size = 3) + # Replace geom_text with geom_text_repel
  #geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2, label = species),
  #                size = 3) + # Replace geom_text with geom_text_repel
  
  #stat_ellipse(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, group = Month), 
  #            linetype = "solid") +
  scale_color_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_shape_manual(labels=c('Control', 'Fill Sand', 'Preemergent','Soil Inversion'),
                     values=c(22, 23, 24, 25)) +
  annotate("text", x = -1, y = 1,
           label = paste0("Stress: ", format(MDS$stress, digits = 2)),
           hjust = 0, size = 8) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="none")+
  labs(x = "MDS1", y = "MDS2", color = "", fill = "", shape = "",
       title = "Post-treatment") 
Post_NMDS_graph

########################## Adonis Analysis #####################################
# Perform Adonis analysis (Permutational MANOVA) to test for significant
# differences in community structure between treatments.
# We use the Bray-Curtis distance matrix and specify the data frame containing the 'Treatment' variable.
# The data frame provided to `adonis2` must have rows in the same order as the distance matrix.
adon.results <- adonis2(Veg_Spp ~ Treatment, data = metadata, method = "bray", permutations = 999)
print("Adonis Results:")
print(adon.results)

# Save the Adonis results to a CSV file.
write.csv(adon.results, "Tables/post_Seedbank_adonis.csv")

# Perform pairwise Adonis to see which specific treatments are different from each other.
# This is only necessary if the main adonis test is significant.
pairwise.adon.results <- pairwise.adonis2(Veg_Spp ~ Treatment, data = metadata,
                                          permutations = 999, method = "bray")
print("Pairwise Adonis Results:")
print(pairwise.adon.results)

# Save the pairwise Adonis results to a CSV file.
write.csv(pairwise.adon.results, "Tables/post_Seedbank_pairwise_adonis.csv")

################################################################################
################################################################################

Pre_Post_Seedbank_NMDS = ggarrange(Pre_NMDS_graph, Post_NMDS_graph, 
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          label.y = 1)
Pre_Post_Seedbank_NMDS


ggsave("Figures/Pre_Post_Seedbank_NMDS.tiff", dpi = 100, width = 10, height = 6)

