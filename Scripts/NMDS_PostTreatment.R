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
                      "devtools", "knitr", "tables", "openxlsx", "labdsv", "permute")
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
library(permute)
library(labdsv)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Load necessary libraries
library(tidyverse)
library(labdsv) # Ensure labdsv is loaded for matrify



################################################################################
######################## Coverage-cleaning helper function #####################
################################################################################

# Preserve every quadrat record used for plot-level averaging.
# Blank, NA, and N/A cells are treated as zero cover because they represent
# surveyed quadrats in which the species or cover class was absent.
clean_coverage_data <- function(df_long, survey_label) {
  
  cleaned <- df_long %>%
    dplyr::mutate(
      Coverage_Original = Coverage,
      Coverage = trimws(as.character(Coverage)),
      Coverage = dplyr::case_when(
        is.na(Coverage) ~ "0",
        Coverage == "" ~ "0",
        toupper(Coverage) %in% c("NA", "N/A", ".") ~ "0",
        TRUE ~ Coverage
      ),
      Coverage = suppressWarnings(as.numeric(Coverage))
    )
  
  # Any remaining NA values came from unrecognized nonnumeric entries.
  invalid_coverage <- cleaned %>%
    dplyr::filter(is.na(Coverage))
  
  if (nrow(invalid_coverage) > 0) {
    
    invalid_file <- file.path(
      "Tables",
      paste0(
        "Invalid_Coverage_Values_",
        gsub("[^A-Za-z0-9]+", "_", survey_label),
        ".csv"
      )
    )
    
    write.csv(
      invalid_coverage,
      invalid_file,
      row.names = FALSE
    )
    
    stop(
      survey_label,
      " contains ", nrow(invalid_coverage),
      " cover values that could not be converted to numbers. See ",
      invalid_file, "."
    )
  }
  
  # Convert cover classes to percent-cover midpoints.
  cleaned %>%
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
        TRUE ~ Coverage
      )
    ) %>%
    dplyr::select(-Coverage_Original)
}

################################################################################
##################### Plot-level PERMANOVA helper functions ####################
################################################################################

# Treatments were applied to whole plots; therefore, quadrats are subsamples.
# This function averages the quadrats within each plot before PERMANOVA.
# Bare ground remains included as a community variable.

run_plot_level_permanova <- function(df_cleaned, survey_label,
                                     overall_csv, pairwise_xlsx,
                                     nperm = 9999,
                                     expected_quadrats = 10) {
  
  # Check for duplicate records of the same species within a plot and quadrat.
  duplicate_check <- df_cleaned %>%
    dplyr::count(
      Block, Plot, Treatment, Species, Quadrat_ID,
      name = "Number_of_Records"
    ) %>%
    dplyr::filter(Number_of_Records > 1)
  
  if (nrow(duplicate_check) > 0) {
    duplicate_file <- file.path(
      "Tables",
      paste0(
        "Duplicate_Quadrat_Records_",
        gsub("[^A-Za-z0-9]+", "_", survey_label),
        ".csv"
      )
    )
    
    write.csv(duplicate_check, duplicate_file, row.names = FALSE)
    
    stop(
      survey_label,
      " contains duplicate species records within one or more quadrats. ",
      "See ", duplicate_file, " before rerunning the PERMANOVA."
    )
  }
  
  # Confirm that every plot-species combination contains all expected quadrats.
  # Zeros must be retained as true species absences. Blank or NA values removed
  # before this function will reduce the quadrat count and trigger this check.
  quadrat_check <- df_cleaned %>%
    dplyr::group_by(Block, Plot, Treatment, Species) %>%
    dplyr::summarise(
      Quadrats_Observed = dplyr::n_distinct(Quadrat_ID),
      Values_Observed = dplyr::n(),
      .groups = "drop"
    )
  
  incomplete_groups <- quadrat_check %>%
    dplyr::filter(
      Quadrats_Observed != expected_quadrats |
        Values_Observed != expected_quadrats
    )
  
  quadrat_check_file <- file.path(
    "Tables",
    paste0(
      "Quadrat_Count_Check_",
      gsub("[^A-Za-z0-9]+", "_", survey_label),
      ".csv"
    )
  )
  
  write.csv(quadrat_check, quadrat_check_file, row.names = FALSE)
  
  if (nrow(incomplete_groups) > 0) {
    incomplete_file <- file.path(
      "Tables",
      paste0(
        "Incomplete_Plot_Species_Groups_",
        gsub("[^A-Za-z0-9]+", "_", survey_label),
        ".csv"
      )
    )
    
    write.csv(incomplete_groups, incomplete_file, row.names = FALSE)
    
    stop(
      survey_label, " contains ", nrow(incomplete_groups),
      " plot-species combinations without exactly ", expected_quadrats,
      " quadrat values. Determine whether blanks or NA values represent true ",
      "species absences or genuinely missing surveys. See ", incomplete_file, "."
    )
  }
  
  # Average all expected quadrats within each plot for each species/cover class.
  # Because every group passed the checks above, the mean includes all quadrats,
  # including explicit zeros for absences.
  plot_level_long <- df_cleaned %>%
    dplyr::group_by(Block, Plot, Treatment, Species) %>%
    dplyr::summarise(
      Mean_Coverage = mean(Coverage),
      Quadrats_Used = dplyr::n_distinct(Quadrat_ID),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Block = factor(Block),
      Plot = factor(Plot),
      Treatment = droplevels(factor(
        Treatment,
        levels = c("Control", "Fill Sand", "Preemergent", "Soil Inversion")
      )),
      Plot_ID = interaction(Block, Plot, drop = TRUE)
    )
  
  # Create one row per plot and one column per species/ground-cover category.
  plot_level_wide <- plot_level_long %>%
    dplyr::select(Plot_ID, Block, Plot, Treatment, Species, Mean_Coverage) %>%
    tidyr::pivot_wider(
      names_from = Species,
      values_from = Mean_Coverage,
      values_fill = 0
    ) %>%
    dplyr::arrange(Block, Plot)
  
  metadata <- plot_level_wide %>%
    dplyr::select(Plot_ID, Block, Plot, Treatment) %>%
    as.data.frame()
  
  community_matrix <- plot_level_wide %>%
    dplyr::select(-Plot_ID, -Block, -Plot, -Treatment) %>%
    as.data.frame()
  
  rownames(community_matrix) <- as.character(metadata$Plot_ID)
  
  # Remove variables absent from every plot.
  community_matrix <- community_matrix[
    , colSums(community_matrix, na.rm = TRUE) > 0, drop = FALSE
  ]
  
  # Stop if plot-level aggregation did not produce the expected 24 plots.
  if (nrow(metadata) != 24) {
    stop(
      survey_label, " produced ", nrow(metadata),
      " plot-level observations; 24 were expected."
    )
  }
  
  # Overall PERMANOVA: test treatment using unrestricted permutations.
  # Block remains in the metadata only for plot identification and is not
  # included as a model covariate or permutation restriction.
  overall_model <- vegan::adonis2(
    community_matrix ~ Treatment,
    data = metadata,
    method = "bray",
    permutations = nperm,
    by = "terms"
  )
  
  print(overall_model)
  
  overall_table <- as.data.frame(overall_model)
  overall_table$Term <- rownames(overall_table)
  overall_table <- overall_table %>%
    dplyr::select(Term, dplyr::everything())
  
  write.csv(
    overall_table,
    overall_csv,
    row.names = FALSE
  )
  
  # Pairwise PERMANOVA comparisons at the plot level.
  treatment_pairs <- combn(
    levels(droplevels(metadata$Treatment)),
    2,
    simplify = FALSE
  )
  
  pairwise_results <- purrr::map_dfr(
    treatment_pairs,
    function(pair) {
      
      keep <- metadata$Treatment %in% pair
      
      pair_metadata <- droplevels(metadata[keep, , drop = FALSE])
      pair_community <- community_matrix[keep, , drop = FALSE]
      
      pair_model <- vegan::adonis2(
        pair_community ~ Treatment,
        data = pair_metadata,
        method = "bray",
        permutations = nperm,
        by = "terms"
      )
      
      pair_table <- as.data.frame(pair_model)
      
      treatment_row <- pair_table[
        "Treatment",
        ,
        drop = FALSE
      ]
      
      residual_row <- pair_table[
        "Residual",
        ,
        drop = FALSE
      ]
      
      total_row <- pair_table[
        "Total",
        ,
        drop = FALSE
      ]
      
      tibble::tibble(
        Survey = survey_label,
        Comparison = paste(pair, collapse = " vs "),
        
        Treatment_Df = treatment_row$Df,
        Residual_Df = residual_row$Df,
        Total_Df = total_row$Df,
        
        Treatment_SumOfSqs = treatment_row$SumOfSqs,
        Residual_SumOfSqs = residual_row$SumOfSqs,
        Total_SumOfSqs = total_row$SumOfSqs,
        
        R2 = treatment_row$R2,
        F = treatment_row$F,
        P_Value = treatment_row$`Pr(>F)`
      )
    }
  ) %>%
    dplyr::mutate(
      P_Adjusted_Holm = p.adjust(P_Value, method = "holm")
    )
  
  print(pairwise_results)
  
  openxlsx::write.xlsx(
    pairwise_results,
    pairwise_xlsx,
    overwrite = TRUE
  )
  
  invisible(list(
    overall_model = overall_model,
    overall_table = overall_table,
    pairwise_table = pairwise_results,
    community_matrix = community_matrix,
    metadata = metadata
  ))
}


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

# Preserve all ten quadrats. Blank and NA/N/A cells are treated as zero cover.
df_cleaned <- clean_coverage_data(
  df_long = df_long,
  survey_label = "6-month"
)

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
  "Tables/NMDS_FitValues_April.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Remove species absent from every quadrat
Spp_filtered <- Spp[
  , colSums(Spp, na.rm = TRUE) > 0,
  drop = FALSE
]

# Perform Indicator Species Analysis
indicator_results <- labdsv::indval(
  Spp_filtered,
  Treat$Treatment
)

# Extract the indicator-value matrix
indval_df <- as.data.frame(indicator_results$indval)

# Preserve species names
indval_df$Species <- rownames(indval_df)

# Move Species to the first column
indval_df <- indval_df %>%
  dplyr::select(
    Species,
    dplyr::everything()
  )

# Add class and significance results
combined_df <- indval_df %>%
  dplyr::mutate(
    Max_Class = indicator_results$maxcls,
    Indicator_Class = indicator_results$indcls,
    P_Value = indicator_results$pval
  )

# Retain significant indicator taxa
significant_species <- combined_df %>%
  dplyr::filter(P_Value < 0.05) %>%
  dplyr::arrange(P_Value) %>%
  dplyr::mutate(
    Rank = dplyr::row_number()
  )

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Tables/indicator_species_analysis_results_April.csv", 
  row.names = FALSE)

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
  file = "Tables/SIMPER_April.csv", 
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

# Perform plot-level PERMANOVA to test treatment effects
PERMANOVA_April <- run_plot_level_permanova(
  df_cleaned = df_cleaned,
  survey_label = "6-month",
  overall_csv = "Tables/adonis_April.csv",
  pairwise_xlsx = "Tables/pairwise_adonis_same_sheet_April.xlsx",
  nperm = 9999
)
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

# Preserve all ten quadrats. Blank and NA/N/A cells are treated as zero cover.
df_cleaned <- clean_coverage_data(
  df_long = df_long,
  survey_label = "3-month"
)

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
  "Tables/NMDS_FitValues_Jan.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Remove species absent from every quadrat
Spp_filtered <- Spp[
  , colSums(Spp, na.rm = TRUE) > 0,
  drop = FALSE
]

# Perform Indicator Species Analysis
indicator_results <- labdsv::indval(
  Spp_filtered,
  Treat$Treatment
)

# Extract the indicator-value matrix
indval_df <- as.data.frame(indicator_results$indval)

# Preserve species names
indval_df$Species <- rownames(indval_df)

# Move Species to the first column
indval_df <- indval_df %>%
  dplyr::select(
    Species,
    dplyr::everything()
  )

# Add class and significance results
combined_df <- indval_df %>%
  dplyr::mutate(
    Max_Class = indicator_results$maxcls,
    Indicator_Class = indicator_results$indcls,
    P_Value = indicator_results$pval
  )

# Retain significant indicator taxa
significant_species <- combined_df %>%
  dplyr::filter(P_Value < 0.05) %>%
  dplyr::arrange(P_Value) %>%
  dplyr::mutate(
    Rank = dplyr::row_number()
  )

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Tables/indicator_species_analysis_results_Jan.csv", 
  row.names = FALSE)

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
  file = "Tables/SIMPER_Jan.csv", 
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

# Perform plot-level PERMANOVA to test treatment effects
PERMANOVA_Jan <- run_plot_level_permanova(
  df_cleaned = df_cleaned,
  survey_label = "3-month",
  overall_csv = "Tables/adonis_Jan.csv",
  pairwise_xlsx = "Tables/pairwise_adonis_same_sheet_Jan.xlsx",
  nperm = 9999
)
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

# Preserve all ten quadrats. Blank and NA/N/A cells are treated as zero cover.
df_cleaned <- clean_coverage_data(
  df_long = df_long,
  survey_label = "1-month"
)

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
  "Tables/NMDS_FitValues_Nov.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Remove species absent from every quadrat
Spp_filtered <- Spp[
  , colSums(Spp, na.rm = TRUE) > 0,
  drop = FALSE
]

# Perform Indicator Species Analysis
indicator_results <- labdsv::indval(
  Spp_filtered,
  Treat$Treatment
)

# Extract the indicator-value matrix
indval_df <- as.data.frame(indicator_results$indval)

# Preserve species names
indval_df$Species <- rownames(indval_df)

# Move Species to the first column
indval_df <- indval_df %>%
  dplyr::select(
    Species,
    dplyr::everything()
  )

# Add class and significance results
combined_df <- indval_df %>%
  dplyr::mutate(
    Max_Class = indicator_results$maxcls,
    Indicator_Class = indicator_results$indcls,
    P_Value = indicator_results$pval
  )

# Retain significant indicator taxa
significant_species <- combined_df %>%
  dplyr::filter(P_Value < 0.05) %>%
  dplyr::arrange(P_Value) %>%
  dplyr::mutate(
    Rank = dplyr::row_number()
  )

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Tables/indicator_species_analysis_results_Nov.csv", 
  row.names = FALSE)

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
  file = "Tables/SIMPER_Nov.csv", 
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

# Perform plot-level PERMANOVA to test treatment effects
PERMANOVA_Nov <- run_plot_level_permanova(
  df_cleaned = df_cleaned,
  survey_label = "1-month",
  overall_csv = "Tables/adonis_Nov.csv",
  pairwise_xlsx = "Tables/pairwise_adonis_same_sheet_Nov.xlsx",
  nperm = 9999
)
################################################################################
################################################################################
NMDS_Combined = ggarrange(Nov_NMDS_graph, Nov_NMDS_graph_Spp,
                          NMDS_graph_Jan , NMDS_graph_Jan_Spp,
                          NMDS_graph_April, NMDS_graph_April_Spp,
                          nrow = 3, ncol = 2, common.legend = TRUE,
                          legend = "bottom")
NMDS_Combined
ggsave("Figures/NMDS_Combined.tiff", dpi = 100, width = 18, height = 24)