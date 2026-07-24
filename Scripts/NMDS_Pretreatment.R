################################################################################
################################################################################
#########################    Site Prep Study      ##############################
#########################    NMDS - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2024 - 2025         ##############################
################################################################################
################################################################################

rm(list=ls(all=TRUE))
cat("\014")

list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont",
                      "ggrepel", "ggsignif", "multcompView", "ggpubr",
                      "rstatix", "rmarkdown", "labdsv", "pairwiseAdonis",
                      "devtools", "knitr", "tables", "openxlsx")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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
library(labdsv)
library(devtools)
library(knitr)
library(tables)
library(openxlsx)

if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)

Data = read.csv("Data/Pretreatment Data.csv", stringsAsFactors = FALSE)

################################################################################
################################### June  ######################################
################################################################################

# Blanks and NA-style entries represent absence and are converted to zero
df_cleaned <- Data %>%
  mutate(
    Coverage = trimws(as.character(Coverage)),
    Coverage = case_when(
      is.na(Coverage) ~ "0",
      Coverage == "" ~ "0",
      toupper(Coverage) %in% c("NA", "N/A", ".") ~ "0",
      TRUE ~ Coverage
    ),
    Coverage = suppressWarnings(as.numeric(Coverage))
  )

if (any(is.na(df_cleaned$Coverage))) {
  stop("Some Coverage values could not be converted to numeric.")
}

# Create the quadrat-level unique ID
df_cleaned <- df_cleaned %>%
  mutate(Unique_ID = paste(Block, Plot, Quadrat, sep = '-'))

# Reclassify coverage data from 1-10 scale to percent cover
df_cleaned <- df_cleaned %>%
  mutate(Coverage = case_when(
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
    TRUE ~ NA_real_
  ))

if (any(is.na(df_cleaned$Coverage))) {
  stop("Coverage values outside the expected 0-10 scale were found.")
}

if (any(is.na(df_cleaned$Treatment))) {
  stop("Treatment names do not match the expected factor levels.")
}

################################################################################
######################## Average Quadrats Within Plots ##########################
################################################################################

# Treatment was applied at the plot level, so quadrats are subsamples.
# First complete absent species as zero, then average the ten quadrats per plot.
Plot_Long <- df_cleaned %>%
  group_by(Block, Plot, Quadrat, Treatment, Species) %>%
  summarise(Coverage = sum(Coverage, na.rm = TRUE), .groups = "drop") %>%
  complete(
    nesting(Block, Plot, Quadrat, Treatment),
    Species,
    fill = list(Coverage = 0)
  ) %>%
  group_by(Block, Plot, Treatment, Species) %>%
  summarise(
    Coverage = mean(Coverage, na.rm = TRUE),
    Number_of_Quadrats = n_distinct(Quadrat),
    .groups = "drop"
  )

if (any(Plot_Long$Number_of_Quadrats != 10)) {
  dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
  write.csv(
    Plot_Long %>% filter(Number_of_Quadrats != 10),
    "Tables/Pretreatment_Incomplete_Plot_Species.csv",
    row.names = FALSE
  )
  stop("One or more plot-species combinations do not contain ten quadrats.")
}

# Create one row per plot and one column per species
Plot_Wide <- Plot_Long %>%
  dplyr::select(Block, Plot, Treatment, Species, Coverage) %>%
  pivot_wider(
    names_from = Species,
    values_from = Coverage,
    values_fill = 0
  ) %>%
  arrange(Treatment, Block, Plot) %>%
  mutate(Unique_ID = paste(Block, Plot, sep = '-'))

Treat <- Plot_Wide %>%
  dplyr::select(Unique_ID, Block, Plot, Treatment)

Spp <- Plot_Wide %>%
  dplyr::select(-Unique_ID, -Block, -Plot, -Treatment) %>%
  as.data.frame()

rownames(Spp) <- Treat$Unique_ID

# Remove species absent from every plot
Spp <- Spp[, colSums(Spp, na.rm = TRUE) > 0, drop = FALSE]

if (any(rowSums(Spp, na.rm = TRUE) == 0)) {
  stop("At least one plot has zero total vegetation cover.")
}

Veg_Spp = vegdist(Spp, method = 'bray')

################################################################################
############################ NMDS Analysis ######################################
################################################################################

NMDS.scree <- function(x) {
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = FALSE, k = 1)$stress),
       xlim = c(1, 10), ylim = c(0, 0.30), xlab = "# of Dimensions",
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = FALSE, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp)

set.seed(123)
MDS = metaMDS(Spp, distance = 'bray', k=3, trymax = 200,
              autotransform = FALSE, trace = FALSE)
MDS$stress
stressplot(MDS)
goodness(MDS)

species.scores <- as.data.frame(wascores(MDS$points, Spp))
species.scores$species <- rownames(species.scores)

NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,
                  Block = Treat$Block, Plot = Treat$Plot,
                  Unique_ID = Treat$Unique_ID)

################################################################################
############################# NMDS Graph ########################################
################################################################################

NMDS_graph = ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
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
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = paste0("Stress: ", round(MDS$stress, 3)),
    hjust = -0.7,
    vjust = 1.3,
    size = 6,
    fontface = "bold"
  ) +
  labs(x = "MDS1", y = "MDS2", color = "Treatment",
       fill = "Treatment", shape = "Treatment")
NMDS_graph

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures/pretreat_NMDS_graph.tiff", plot = NMDS_graph,
       dpi = 300, width = 8, height = 6, compression = "lzw")

################################################################################
############################ PERMANOVA ##########################################
################################################################################

# Plot-level PERMANOVA; block is included because the experiment was blocked
set.seed(123)
adon.results <- adonis2(
  Spp ~ Block + Treatment,
  data = Treat,
  method = "bray",
  permutations = 9999,
  by = "terms"
)
print(adon.results)

dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(adon.results),
          "Tables/Pretreat_adonis.csv", row.names = TRUE)

################################################################################
######################## Multivariate Dispersion ###############################
################################################################################

dispersion_model <- betadisper(Veg_Spp, Treat$Treatment, type = "centroid")
set.seed(123)
dispersion_results <- permutest(dispersion_model, permutations = 9999)
print(dispersion_results)

write.csv(as.data.frame(dispersion_results$tab),
          "Tables/Pretreat_dispersion.csv", row.names = TRUE)

################################################################################
######################## Pairwise PERMANOVA ####################################
################################################################################

# Retained for consistency, but do not interpret these if the overall
# pretreatment Treatment effect is nonsignificant.
pairwise.adonis <- pairwise.adonis2(
  Spp ~ Treatment,
  data = Treat,
  permutations = 9999
)
pairwise.adonis

wb <- createWorkbook()
addWorksheet(wb, "All_Tables")
start_row <- 1

for (name in names(pairwise.adonis)) {
  writeData(wb, sheet = "All_Tables", x = name,
            startRow = start_row, colNames = FALSE)
  start_row <- start_row + 1
  
  if (is.data.frame(pairwise.adonis[[name]])) {
    writeData(wb, sheet = "All_Tables",
              x = pairwise.adonis[[name]], startRow = start_row)
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    writeData(wb, sheet = "All_Tables",
              x = pairwise.adonis[[name]],
              startRow = start_row, colNames = FALSE)
    start_row <- start_row + 2
  }
}

saveWorkbook(wb, "Tables/Pretreat_pairwise_adonis_same_sheet.xlsx",
             overwrite = TRUE)