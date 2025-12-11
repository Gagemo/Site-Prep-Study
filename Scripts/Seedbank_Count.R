# Load necessary packages
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)
library(emmeans)
library(multcompView)
library(multcomp)

# Set the file name for the uploaded data
# The file "Seedbank.csv" has been made available
file_path <- "Data/Seedbank.csv"

# Read in the Seedbank data
# The file has a different column structure from your previous script, so we'll read it directly.
Data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Filter the data for "Post-Treatment" month only, as this is when the effects of the treatments would be visible.
post_treatment_data <- Data %>%
  filter(Month == "Post-Treatment")

# The species columns start from column 6 to the end.
# We'll calculate the total number of germinated seeds for each sample.
post_germination_data <- post_treatment_data %>%
  rowwise() %>%
  # Sum all species counts for each row (starting from the 6th column)
  mutate(TotalGermination = sum(c_across(Ambrosia.artemisiifolia:Triodanis.perfoliata))) %>%
  ungroup()

# Ensure Treatment is a factor with the correct levels for analysis
post_germination_data$Treatment <- factor(post_germination_data$Treatment,
                                          levels = c("Control",
                                                     "Fill Sand",
                                                     "Preemergent",
                                                     "Soil Inversion"))

# Perform an ANOVA using the rstatix package for better compatibility with ggpubr
post_anova <- post_germination_data %>%
  anova_test(TotalGermination ~ Treatment)

# Display ANOVA results in the console
cat("--------------------------------------------------\n")
cat("ANOVA Table for Post-Treatment Seed Germination by Treatment\n")
cat("--------------------------------------------------\n")
print(post_anova)

# Add a compact letter display to the plot
# Run a post-hoc analysis to get the groups
post_tukey_cld <- emmeans(aov(TotalGermination ~ Treatment, data = post_germination_data), ~Treatment) %>%
  cld(
    adjust = "Tukey",
    Letters = letters,
    decreasing = TRUE,
    alpha = 0.05
  ) %>%
  # Add the letters to the original data frame
  mutate(label = .group)

# Get the maximum y-value for each treatment group for positioning the labels
post_max_y_values <- post_germination_data %>%
  group_by(Treatment) %>%
  summarise(max_y = max(TotalGermination))

# Join the max y-values to the tukey_cld data frame
post_tukey_cld <- left_join(post_tukey_cld, post_max_y_values, by = "Treatment")

# Create a boxplot to visualize the results
# We'll use your previous color palette for consistency
treat_colors <- c("Control" = "yellow",
                  "Fill Sand" = "orange",
                  "Preemergent" = "#66CC00",
                  "Soil Inversion" = "#CC66CC")

# Generate the plot
post <- ggplot(post_germination_data, aes(x = Treatment, y = TotalGermination, fill = Treatment)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(size = 3, alpha = 0.7, width = 0.2, show.legend = FALSE) +
  # Add the compact letters
  geom_text(data = post_tukey_cld, aes(x = Treatment, y = max_y + 2, label = label),
            inherit.aes = FALSE, vjust = 0, hjust = 0.5, size = 8,
            fontface = "bold") +
  # Use the anova_results object to create the subtitle
  labs(subtitle = get_test_label(post_anova, detailed = TRUE)) +
  ylim(0, 250) +
  scale_fill_manual(values = treat_colors) +
  scale_color_manual(values = treat_colors) +
  labs(x = "", y = "Total germination count",
       title = "Post-treatment") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
        text = element_text(size = 20),
        axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 20, face = "bold", color = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y =element_blank(),
        legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.position = "bottom")

# Filter the data for "Pre-Treatment" month
pre_treatment_data <- Data %>%
  filter(Month == "Pre-Treatment")

# Calculate the total number of germinated seeds for each sample.
pre_germination_data <- pre_treatment_data %>%
  rowwise() %>%
  mutate(TotalGermination = sum(c_across(Ambrosia.artemisiifolia:Triodanis.perfoliata))) %>%
  ungroup()

# Ensure Treatment is a factor with the correct levels for analysis
pre_germination_data$Treatment <- factor(pre_germination_data$Treatment,
                                         levels = c("Control",
                                                    "Fill Sand",
                                                    "Preemergent",
                                                    "Soil Inversion"))

# Perform an ANOVA using the rstatix package
pre_anova <- pre_germination_data %>%
  anova_test(TotalGermination ~ Treatment)

# Display ANOVA results in the console
cat("--------------------------------------------------\n")
cat("ANOVA Table for Pre-Treatment Seed Germination by Treatment\n")
cat("--------------------------------------------------\n")
print(pre_anova)

# Add a compact letter display to the plot
pre_tukey_cld <- emmeans(aov(TotalGermination ~ Treatment, data = pre_germination_data), ~Treatment) %>%
  cld(
    adjust = "Tukey",
    Letters = letters,
    decreasing = TRUE,
    alpha = 0.05
  ) %>%
  mutate(label = .group)

# Get the maximum y-value for each treatment group for positioning the labels
pre_max_y_values <- pre_germination_data %>%
  group_by(Treatment) %>%
  summarise(max_y = max(TotalGermination))

# Join the max y-values to the tukey_cld data frame
pre_tukey_cld <- left_join(pre_tukey_cld, pre_max_y_values, by = "Treatment")

# Create a boxplot for pre-treatment
pre <- ggplot(pre_germination_data, aes(x = Treatment, y = TotalGermination, fill = Treatment)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(size = 3, alpha = 0.7, width = 0.2, show.legend = FALSE) +
  geom_text(data = pre_tukey_cld, aes(x = Treatment, y = max_y + 2, label = label),
            inherit.aes = FALSE, vjust = 0, hjust = 0.5, size = 8,
            fontface = "bold") +
  # Use the anova_results object to create the subtitle
  labs(subtitle = get_test_label(pre_anova, detailed = TRUE)) +
  scale_fill_manual(values = treat_colors) +
  scale_color_manual(values = treat_colors) +
  ylim(0, 250) +
  labs(x = "", y = "Total germination count",
       title = "Pre-treatment") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
        text = element_text(size = 20),
        axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 20, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 20, face = "bold", color = "black"),
        axis.text.y = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.position = "bottom")

# Arrange the two plots side-by-side
seedbank_germ_count <- ggarrange(
  pre, 
  post, 
  nrow = 1, 
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)

# Add the post-hoc test label to the bottom of the combined figure
final_figure <- annotate_figure(
  seedbank_germ_count,
  bottom = text_grob("pwc: Dunn test; p.adjust: Holm",
                     hjust = -2,
                     vjust = -1.8,
                     size = 15)
)

# Print the final arranged figure
print(final_figure)

# Save the figure to a file
ggsave("Figures/seedbank_germ_count.png", plot = final_figure, dpi = 600, 
       width = 18, height = 8)

