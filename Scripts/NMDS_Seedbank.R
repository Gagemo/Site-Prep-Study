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
rm(list = ls(all = TRUE))
cat("\014")

######################### Installs & Loads Packages #############################
list.of.packages <- c(
  "tidyverse", "vegan", "ggrepel", "ggpubr", "openxlsx"
)

new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages)) {
  install.packages(new.packages)
}

invisible(lapply(list.of.packages, library, character.only = TRUE))

########################## User Settings ########################################
data_file <- "Data/Seedbank.csv"
output_table_dir <- "Tables"
output_figure_dir <- "Figures"

nperm <- 9999
expected_samples_per_plot <- 4
set.seed(123)

dir.create(output_table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_figure_dir, showWarnings = FALSE, recursive = TRUE)

########################## Read and Prepare Data ################################
Data <- read.csv(
  data_file,
  check.names = TRUE,
  stringsAsFactors = FALSE
)

required_metadata <- c(
  "Month", "Block", "Plot", "Sample", "Treatment"
)

missing_metadata <- setdiff(required_metadata, names(Data))

if (length(missing_metadata) > 0) {
  stop(
    "The seedbank dataset is missing required columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

Data$Treatment <- factor(
  Data$Treatment,
  levels = c(
    "Control",
    "Fill Sand",
    "Preemergent",
    "Soil Inversion"
  )
)

Data$Month <- factor(
  Data$Month,
  levels = c(
    "Pre-Treatment",
    "Post-Treatment"
  )
)

if (any(is.na(Data$Treatment))) {
  stop(
    "One or more Treatment values do not match the expected treatment names."
  )
}

if (any(is.na(Data$Month))) {
  stop(
    "One or more Month values do not match 'Pre-Treatment' or 'Post-Treatment'."
  )
}

# Identify species columns as every column that is not plot/sample metadata.
species_columns <- setdiff(
  names(Data),
  c("Month", "Block", "Plot", "Sample", "Treatment", "Unique_ID")
)

if (length(species_columns) == 0) {
  stop("No species columns were detected in the seedbank dataset.")
}

# Convert all species values to numeric and treat blank/NA/N/A/'.' as zero.
Data <- Data %>%
  mutate(
    across(
      all_of(species_columns),
      ~ {
        x <- trimws(as.character(.x))
        x[
          is.na(x) |
            x == "" |
            toupper(x) %in% c("NA", "N/A", ".")
        ] <- "0"
        suppressWarnings(as.numeric(x))
      }
    )
  )

invalid_species_values <- Data %>%
  mutate(Row_Number = row_number()) %>%
  pivot_longer(
    cols = all_of(species_columns),
    names_to = "Species",
    values_to = "Value"
  ) %>%
  filter(is.na(Value))

if (nrow(invalid_species_values) > 0) {
  invalid_file <- file.path(
    output_table_dir,
    "Seedbank_Invalid_Species_Values.csv"
  )
  
  write.csv(
    invalid_species_values,
    invalid_file,
    row.names = FALSE
  )
  
  stop(
    "One or more seedbank values could not be converted to numeric. ",
    "See ", invalid_file, "."
  )
}

# Confirm that no abundance values are negative.
negative_values <- Data %>%
  pivot_longer(
    cols = all_of(species_columns),
    names_to = "Species",
    values_to = "Value"
  ) %>%
  filter(Value < 0)

if (nrow(negative_values) > 0) {
  negative_file <- file.path(
    output_table_dir,
    "Seedbank_Negative_Values.csv"
  )
  
  write.csv(
    negative_values,
    negative_file,
    row.names = FALSE
  )
  
  stop(
    "Negative seedbank abundance values were found. See ",
    negative_file,
    "."
  )
}

########################## Plot-Level Aggregation ###############################
# Treatments were applied to plots. The four soil seedbank samples collected
# within each plot are subsamples and are averaged before treatment comparisons.

duplicate_samples <- Data %>%
  count(
    Month,
    Block,
    Plot,
    Sample,
    Treatment,
    name = "Number_of_Rows"
  ) %>%
  filter(Number_of_Rows > 1)

if (nrow(duplicate_samples) > 0) {
  duplicate_file <- file.path(
    output_table_dir,
    "Seedbank_Duplicate_Samples.csv"
  )
  
  write.csv(
    duplicate_samples,
    duplicate_file,
    row.names = FALSE
  )
  
  stop(
    "Duplicate seedbank sample records were found. See ",
    duplicate_file,
    "."
  )
}

sample_count_check <- Data %>%
  group_by(
    Month,
    Block,
    Plot,
    Treatment
  ) %>%
  summarise(
    Number_of_Samples = n_distinct(Sample),
    Number_of_Rows = n(),
    .groups = "drop"
  )

incomplete_plots <- sample_count_check %>%
  filter(
    Number_of_Samples != expected_samples_per_plot |
      Number_of_Rows != expected_samples_per_plot
  )

if (nrow(incomplete_plots) > 0) {
  incomplete_file <- file.path(
    output_table_dir,
    "Seedbank_Incomplete_Plots.csv"
  )
  
  write.csv(
    incomplete_plots,
    incomplete_file,
    row.names = FALSE
  )
  
  stop(
    "One or more plot-survey combinations do not contain exactly ",
    expected_samples_per_plot,
    " seedbank samples. See ",
    incomplete_file,
    "."
  )
}

Plot_Data <- Data %>%
  group_by(
    Month,
    Block,
    Plot,
    Treatment
  ) %>%
  summarise(
    across(
      all_of(species_columns),
      ~ mean(.x, na.rm = TRUE)
    ),
    Number_of_Samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  mutate(
    Plot_ID = interaction(Block, Plot, drop = TRUE),
    Unique_ID = paste(Month, Block, Plot, sep = "-")
  )

# Confirm expected plot replication.
plot_replication <- Plot_Data %>%
  count(
    Month,
    Treatment,
    name = "Number_of_Plots"
  )

write.csv(
  plot_replication,
  file.path(output_table_dir, "Seedbank_Plot_Replication.csv"),
  row.names = FALSE
)

if (any(plot_replication$Number_of_Plots != 6)) {
  warning(
    "At least one treatment-by-survey group does not contain six plots. ",
    "Review Tables/Seedbank_Plot_Replication.csv."
  )
}

write.csv(
  Plot_Data,
  file.path(output_table_dir, "Seedbank_Plot_Level_Data.csv"),
  row.names = FALSE
)

########################## Analysis Helper Functions ############################
analyze_seedbank_survey <- function(
    plot_data,
    survey_name,
    file_stub,
    nperm = 9999
) {
  
  survey_data <- plot_data %>%
    filter(Month == survey_name) %>%
    arrange(Treatment, Block, Plot)
  
  metadata <- survey_data %>%
    dplyr::select(
      Unique_ID,
      Month,
      Block,
      Plot,
      Plot_ID,
      Treatment,
      Number_of_Samples
    )
  
  species_matrix <- survey_data %>%
    dplyr::select(all_of(species_columns)) %>%
    as.data.frame()
  
  rownames(species_matrix) <- metadata$Unique_ID
  
  # Remove species absent from every plot during this survey.
  all_zero_species <- names(species_matrix)[
    colSums(species_matrix, na.rm = TRUE) == 0
  ]
  
  if (length(all_zero_species) > 0) {
    write.csv(
      data.frame(Species = all_zero_species),
      file.path(
        output_table_dir,
        paste0(file_stub, "_All_Zero_Species.csv")
      ),
      row.names = FALSE
    )
    
    species_matrix <- species_matrix[
      ,
      colSums(species_matrix, na.rm = TRUE) > 0,
      drop = FALSE
    ]
  }
  
  if (ncol(species_matrix) < 2) {
    stop(
      survey_name,
      " has fewer than two nonzero species after filtering."
    )
  }
  
  if (any(rowSums(species_matrix, na.rm = TRUE) == 0)) {
    zero_plot_file <- file.path(
      output_table_dir,
      paste0(file_stub, "_Zero_Total_Plots.csv")
    )
    
    write.csv(
      metadata[
        rowSums(species_matrix, na.rm = TRUE) == 0,
        ,
        drop = FALSE
      ],
      zero_plot_file,
      row.names = FALSE
    )
    
    stop(
      survey_name,
      " contains one or more plots with zero total emergence. See ",
      zero_plot_file,
      "."
    )
  }
  
  ########################## Bray-Curtis and NMDS ###############################
  bray_distance <- vegan::vegdist(
    species_matrix,
    method = "bray"
  )
  
  set.seed(123)
  MDS <- vegan::metaMDS(
    species_matrix,
    distance = "bray",
    k = 3,
    trymax = 200,
    autotransform = FALSE,
    trace = FALSE
  )
  
  site_scores <- as.data.frame(
    vegan::scores(MDS, display = "sites")
  ) %>%
    rownames_to_column("Unique_ID") %>%
    left_join(metadata, by = "Unique_ID")
  
  species_scores <- as.data.frame(
    vegan::scores(MDS, display = "species")
  ) %>%
    rownames_to_column("Species")
  
  ########################## Overall PERMANOVA ##################################
  set.seed(123)
  overall_model <- vegan::adonis2(
    species_matrix ~ Treatment,
    data = metadata,
    method = "bray",
    permutations = nperm,
    by = "terms"
  )
  
  overall_table <- as.data.frame(overall_model) %>%
    rownames_to_column("Term") %>%
    mutate(
      Survey = survey_name,
      .before = 1
    )
  
  write.csv(
    overall_table,
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_PERMANOVA.csv")
    ),
    row.names = FALSE
  )
  
  ########################## Pairwise PERMANOVA #################################
  treatment_pairs <- combn(
    levels(droplevels(metadata$Treatment)),
    2,
    simplify = FALSE
  )
  
  pairwise_results <- purrr::map_dfr(
    treatment_pairs,
    function(pair) {
      
      pair_rows <- metadata$Treatment %in% pair
      
      pair_metadata <- metadata[
        pair_rows,
        ,
        drop = FALSE
      ] %>%
        mutate(
          Treatment = droplevels(Treatment)
        )
      
      pair_species <- species_matrix[
        pair_rows,
        ,
        drop = FALSE
      ]
      
      stopifnot(
        nrow(pair_metadata) == 12,
        nlevels(pair_metadata$Treatment) == 2
      )
      
      set.seed(123)
      pair_model <- vegan::adonis2(
        pair_species ~ Treatment,
        data = pair_metadata,
        method = "bray",
        permutations = nperm,
        by = "terms"
      )
      
      pair_table <- as.data.frame(pair_model)
      
      if (!all(c("Treatment", "Residual", "Total") %in% rownames(pair_table))) {
        stop(
          "Expected Treatment, Residual, and Total rows were not found for ",
          paste(pair, collapse = " vs "),
          "."
        )
      }
      
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
      
      tibble(
        Survey = survey_name,
        Comparison = paste(pair, collapse = " vs "),
        Treatment_Df = treatment_row$Df,
        Residual_Df = residual_row$Df,
        Total_Df = total_row$Df,
        SumOfSqs = treatment_row$SumOfSqs,
        R2 = treatment_row$R2,
        F = treatment_row$F,
        P_Value = treatment_row$`Pr(>F)`
      )
    }
  ) %>%
    mutate(
      Holm_Adjusted_P = p.adjust(
        P_Value,
        method = "holm"
      )
    )
  
  write.csv(
    pairwise_results,
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_Pairwise_PERMANOVA.csv")
    ),
    row.names = FALSE
  )
  
  openxlsx::write.xlsx(
    list(
      Overall_PERMANOVA = overall_table,
      Pairwise_PERMANOVA = pairwise_results
    ),
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_PERMANOVA_Results.xlsx")
    ),
    overwrite = TRUE
  )
  
  ########################## Multivariate Dispersion ############################
  dispersion_model <- vegan::betadisper(
    bray_distance,
    metadata$Treatment,
    type = "centroid"
  )
  
  dispersion_anova <- as.data.frame(
    anova(dispersion_model)
  ) %>%
    rownames_to_column("Term") %>%
    mutate(
      Survey = survey_name,
      .before = 1
    )
  
  set.seed(123)
  dispersion_permutation <- vegan::permutest(
    dispersion_model,
    permutations = nperm
  )
  
  dispersion_permutation_table <- as.data.frame(
    dispersion_permutation$tab
  ) %>%
    rownames_to_column("Term") %>%
    mutate(
      Survey = survey_name,
      .before = 1
    )
  
  dispersion_distances <- tibble(
    Survey = survey_name,
    Unique_ID = metadata$Unique_ID,
    Treatment = metadata$Treatment,
    Distance_to_Centroid = dispersion_model$distances
  )
  
  dispersion_summary <- dispersion_distances %>%
    group_by(
      Survey,
      Treatment
    ) %>%
    summarise(
      Mean_Distance = mean(Distance_to_Centroid),
      SD_Distance = sd(Distance_to_Centroid),
      SE_Distance = SD_Distance / sqrt(n()),
      Number_of_Plots = n(),
      .groups = "drop"
    )
  
  write.csv(
    dispersion_anova,
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_Dispersion_ANOVA.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    dispersion_permutation_table,
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_Dispersion_Permutation.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    dispersion_summary,
    file.path(
      output_table_dir,
      paste0(file_stub, "_Seedbank_Dispersion_Summary.csv")
    ),
    row.names = FALSE
  )
  
  ########################## NMDS Plot ##########################################
  treatment_colors <- c(
    "Control" = "yellow",
    "Fill Sand" = "orange",
    "Preemergent" = "#66CC00",
    "Soil Inversion" = "#CC66CC"
  )
  
  nmds_plot <- ggplot(
    site_scores,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = Treatment,
      shape = Treatment
    )
  ) +
    geom_point(
      alpha = 0.8,
      size = 5
    ) +
    scale_color_manual(values = treatment_colors) +
    scale_fill_manual(values = treatment_colors) +
    scale_shape_manual(
      values = c(22, 23, 24, 25)
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = paste0(
        "Stress: ",
        format(MDS$stress, digits = 3)
      ),
      hjust = -0.05,
      vjust = 1.2,
      size = 5
    ) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(
        color = "black",
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(
        size = 18,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 18,
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      legend.title = element_blank()
    ) +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      title = survey_name
    )
  
  return(
    list(
      Survey = survey_name,
      Metadata = metadata,
      Species_Matrix = species_matrix,
      Bray_Distance = bray_distance,
      NMDS_Model = MDS,
      Site_Scores = site_scores,
      Species_Scores = species_scores,
      NMDS_Plot = nmds_plot,
      Overall_PERMANOVA = overall_table,
      Pairwise_PERMANOVA = pairwise_results,
      Dispersion_ANOVA = dispersion_anova,
      Dispersion_Permutation = dispersion_permutation_table,
      Dispersion_Summary = dispersion_summary
    )
  )
}

########################## Run Pretreatment Analysis ############################
Pre_Results <- analyze_seedbank_survey(
  plot_data = Plot_Data,
  survey_name = "Pre-Treatment",
  file_stub = "Pretreatment",
  nperm = nperm
)

Pre_NMDS_graph <- Pre_Results$NMDS_Plot +
  theme(
    legend.position = "none"
  )

print(Pre_Results$Overall_PERMANOVA)
print(Pre_Results$Pairwise_PERMANOVA)
print(Pre_Results$Dispersion_Permutation)
print(Pre_NMDS_graph)

########################## Run Post-treatment Analysis ##########################
Post_Results <- analyze_seedbank_survey(
  plot_data = Plot_Data,
  survey_name = "Post-Treatment",
  file_stub = "Posttreatment",
  nperm = nperm
)

Post_NMDS_graph <- Post_Results$NMDS_Plot +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

print(Post_Results$Overall_PERMANOVA)
print(Post_Results$Pairwise_PERMANOVA)
print(Post_Results$Dispersion_Permutation)
print(Post_NMDS_graph)

########################## Combined Results Files ###############################
Combined_Overall_PERMANOVA <- bind_rows(
  Pre_Results$Overall_PERMANOVA,
  Post_Results$Overall_PERMANOVA
)

Combined_Pairwise_PERMANOVA <- bind_rows(
  Pre_Results$Pairwise_PERMANOVA,
  Post_Results$Pairwise_PERMANOVA
)

Combined_Dispersion <- bind_rows(
  Pre_Results$Dispersion_Summary,
  Post_Results$Dispersion_Summary
)

openxlsx::write.xlsx(
  list(
    Overall_PERMANOVA = Combined_Overall_PERMANOVA,
    Pairwise_PERMANOVA = Combined_Pairwise_PERMANOVA,
    Dispersion_Summary = Combined_Dispersion,
    Pretreatment_Dispersion_Test =
      Pre_Results$Dispersion_Permutation,
    Posttreatment_Dispersion_Test =
      Post_Results$Dispersion_Permutation
  ),
  file.path(
    output_table_dir,
    "Seedbank_Combined_Results.xlsx"
  ),
  overwrite = TRUE
)

########################## Combined NMDS Figure #################################
Pre_Post_Seedbank_NMDS <- ggpubr::ggarrange(
  Pre_NMDS_graph,
  Post_NMDS_graph,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  labels = c("A", "B")
)

print(Pre_Post_Seedbank_NMDS)

ggsave(
  filename = file.path(
    output_figure_dir,
    "Pre_Post_Seedbank_NMDS.tiff"
  ),
  plot = Pre_Post_Seedbank_NMDS,
  dpi = 300,
  width = 12,
  height = 7,
  units = "in",
  compression = "lzw"
)

################################################################################
########################## End of Seedbank Analysis #############################
################################################################################