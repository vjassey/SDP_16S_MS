#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 2B - Determinstic versus stochastic processes
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# This script is adapted from Riddley, M., Hepp, S., Hardeep, F. et al. 
#Differential roles of deterministic and stochastic processes in structuring soil bacterial ecotypes across terrestrial ecosystems. 
#Nat Commun 16, 2337 (2025). https://doi.org/10.1038/s41467-025-57526-x

#See Riddley et al. for details on how to get the input files

# Load required libraries
library(tidyverse)
library(readr)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

load("SDPs_16S.RData")

# Refine format -----------------------------------------------------------

# Define the reformat_output function
reformat_output <- function(df) {
  # Drop the "sig.index" column
  df <- df %>% select(-`sig.index`)
  
  # Replace "index." with "" and dots with "-" in column names
  colnames(df) <- colnames(df) %>%
    str_replace("index\\.", "") %>%
    str_replace("\\.", "-")
  
  # Set index using column names
  rownames(df) <- colnames(df)
  
  return(df)
}

# Read the input CSV file
df_all <- read.delim("BNTI_all.csv", sep=',') #input file

# Reformat the data frame
df_all_refor <- reformat_output(df_all)

# Save the reformatted data frame to a CSV file
write_csv(df_all_refor, "BNTI_all_refor.csv") #input file


# NTI Distribution --------------------------------------------------------

# Read data
df_all <- read.delim("NTI_abun_rownames.csv", sep=',') #input file

# Merge all ecotypes
merged_df <- df_all # %>%
#left_join(df_abun, by = "Sample ID") %>%
#left_join(df_gen, by = "Sample ID") %>%
#left_join(df_rare, by = "Sample ID") %>%
#left_join(df_spe, by = "Sample ID")

colnames(merged_df) <- c("Sample.ID", "All")

# Rename columns
#merged_df <- merged_df %>%
#  rename(All = NTI_all, Abundant = NTI_abun, Generalist = NTI_gen, Rare = NTI_rare, Specialist = NTI_spe) %>%
#   select(Sample.ID = Sample.ID, All, Abundant, Rare, Generalist, Specialist)

# Group NTI values into 3 bins: <-2, >= -2 or <= 2, > 2
#NTI_groups <- list()
#for (col in c("All", "Abundant", "Rare", "Generalist", "Specialist")) {
#  less2 <- merged_df$Sample.ID[merged_df[[col]] < -2]
#  btwn2 <- merged_df$Sample.ID[merged_df[[col]] >= -2 & merged_df[[col]] <= 2]
#  more2 <- merged_df$Sample.ID[merged_df[[col]] > 2]
#  NTI_groups[[col]] <- list(less2 = less2, btwn2 = btwn2, more2 = more2)
#}

# Initialiser la liste NTI_groups
NTI_groups <- list()

# Traiter uniquement la colonne "All"
col <- "All"
less2 <- merged_df$Sample.ID[merged_df[[col]] < -2]
btwn2 <- merged_df$Sample.ID[which(merged_df[[col]] >= -2 & merged_df[[col]] <= 2)]
more2 <- merged_df$Sample.ID[which(merged_df[[col]] > 2)]

# Stocker les résultats dans NTI_groups
NTI_groups[[col]] <- list(
  less2 = as.vector(less2),
  btwn2 = as.vector(btwn2),
  more2 = as.vector(more2)
)


# Count NTI bins
NTI_group_count <- data.frame(
  "less2" = sapply(NTI_groups, function(x) length(x$less2)),
  "btwn2" = sapply(NTI_groups, function(x) length(x$btwn2)),
  "more2" = sapply(NTI_groups, function(x) length(x$more2)),
  row.names = names(NTI_groups)
)

# Calculate percentages
NTI_percent <- NTI_group_count %>%
  mutate(across(everything(), ~ .x / rowSums(NTI_group_count) * 100))

# Plot NTI bins
ggplot(NTI_percent, aes(x = rownames(NTI_percent), y = `less2`, fill = "less2")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = `btwn2`, fill = "btwn2"), stat = "identity") +
  geom_bar(aes(y = `more2`, fill = "more2"), stat = "identity") +
  scale_fill_manual(values = c("less2" = "red", "btwn2" = "blue", "more2" = "green")) +
  labs(title = "NTI Distribution by Ecotypes", y = "Percentage (%)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("../output/NTI_ecotypes.pdf", width = 6, height = 6, dpi = 800)

# Ecosystems ----------------------------------------------------------------

# Read data
eco <- SDPs_16s@sam_data
eco <- eco %>% dplyr::select(Sample, koppen_clim)
colnames(eco) <- c("Sample.ID", "koppen")

# Merge with NTI data
merged_eco <- merged_df %>%
  left_join(eco, by = "Sample.ID")

# Categorize NTI values
merged_eco <- merged_eco %>%
  mutate(Category = case_when(
    All < -2 ~ "NTI < -2",
    All >= -2 & All <= 2 ~ "-2 <= NTI <= 2",
    All > 2 ~ "NTI > 2"
  ))

# Calculate percentage of each category within each ecosystem
category_percentage <- merged_eco %>%
  group_by(koppen, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(koppen) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() %>%
  pivot_wider(names_from = Category, values_from = Percentage) %>%
  arrange(`-2 <= NTI <= 2`)

# Plot NTI distribution by ecosystems
ggplot(category_percentage, aes(x = koppen, y = `-2 <= NTI <= 2`, fill = "-2 <= NTI <= 2")) +
  geom_bar(stat = "identity") +
  #geom_bar(aes(y = `-2 <= NTI <= 2`, fill = "-2 <= NTI <= 2"), stat = "identity") +
  geom_bar(aes(y = `NTI > 2`, fill = "NTI > 2"), stat = "identity") +
  scale_fill_manual(values = c("NTI < -2" = "red", "-2 <= NTI <= 2" = "blue", "NTI > 2" = "green")) +
  labs(title = "NTI Distribution by Ecosystems", y = "Percentage (%)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("../output/NTI_ecosystems.pdf", width = 6, height = 6, dpi = 800)

# BNTI Distribution --------------------------------------------------------

# Preprocess BNTI files
preprocess_file <- function(df_taxa, taxonomy) {
  df_taxa <- df_taxa %>%
    mutate(Samples = rownames(.)) %>%
    pivot_longer(cols = -Samples, names_to = "Sample_ID_y", values_to = taxonomy) %>%
    filter(Samples != Sample_ID_y) %>%
    arrange(Samples, Sample_ID_y) %>%
    rename(Sample_ID_x = Samples)
  return(df_taxa)
}

# Load and melt BNTI data
df_all <- df_all_refor %>% preprocess_file("all_value")

# Stratify BNTI by ecosystem

filter_by_ecosystem <- function(ecosystem) {
  samples <- eco %>% filter(koppen == ecosystem) %>% pull(Sample.ID)
  df_all %>% filter(Sample_ID_x %in% samples & Sample_ID_y %in% samples)
}

df_Temperate <- filter_by_ecosystem("Temperate")
df_CCT <- filter_by_ecosystem("Continental cold/Tundra")
df_CHS <- filter_by_ecosystem("Continental hot summer")


# Stochastic Processes -----------------------------------------------------

stochastic_process <- function(df_all, taxonomy_value, df_rc_bray) {
  # Filter the data for stochastic processes
  df_stochastic <- df_all %>%
    filter(.data[[taxonomy_value]] >= -2 & .data[[taxonomy_value]] <= 2)
  
  # Remove "Unnamed: 0" column if it exists
  if ("Unnamed: 0" %in% colnames(df_rc_bray)) {
    df_rc_bray <- df_rc_bray %>% select(-`Unnamed: 0`)
  }
  
  # Merge the data frames
  rc_bray <- df_rc_bray %>%
    inner_join(df_stochastic, by = c("Sample_ID_x", "Sample_ID_y"))
  
  return(rc_bray)
}

# Load RCBray results
df_bray_curtis_all <- read.delim("RCbray_Results_all.csv", sep=",")
colnames(df_bray_curtis_all) <- c("ID", "Sample_ID_x", "Sample_ID_y", "RCbray")

# Apply stochastic process
rc_bray_all <- stochastic_process(df_all, 'all_value', df_bray_curtis_all) 


# Filter RCBray by ecosystem
df_bray_curtis_TEMP <- df_bray_curtis_all %>% filter(Sample_ID_x %in% df_Temperate$Sample_ID_x & Sample_ID_y %in% df_Temperate$Sample_ID_y)
df_bray_curtis_CCT <- df_bray_curtis_all %>% filter(Sample_ID_x %in% df_CCT$Sample_ID_x & Sample_ID_y %in% df_CCT$Sample_ID_y)
df_bray_curtis_CHS <- df_bray_curtis_all %>% filter(Sample_ID_x %in% df_CHS$Sample_ID_x & Sample_ID_y %in% df_CHS$Sample_ID_y)

rc_bray_TEMP <- stochastic_process(df_Temperate, 'all_value', df_bray_curtis_TEMP)
rc_bray_CCT <- stochastic_process(df_CCT, 'all_value', df_bray_curtis_CCT)
rc_bray_CHS <- stochastic_process(df_CHS, 'all_value', df_bray_curtis_CHS)


# Deterministic Processes --------------------------------------------------

deterministic_processes <- function(df_taxa, taxonomy_value) {
  df_deterministic <- df_taxa
  
  df_deterministic_hetero <- df_deterministic %>% filter(.data[[taxonomy_value]] > 2)
  df_deterministic_homo <- df_deterministic %>% filter(.data[[taxonomy_value]] < -2)
  
  return(list(df_deterministic_hetero, df_deterministic_homo))
}

df_all_deterministic <- deterministic_processes(df_all, 'all_value')

df_TEMP_deterministic <- deterministic_processes(df_Temperate, 'all_value')
df_CCT_deterministic <- deterministic_processes(df_CCT, 'all_value')
df_CHS_deterministic <- deterministic_processes(df_CHS, 'all_value')


# Pie Chart for All Taxa, Ecotypes, and Ecosystems --------------------------

library(ggplot2)
library(dplyr)
library(patchwork)

# Sizes (as you already computed)
det_hetero  <- nrow(df_CHS_deterministic[[1]])
det_homo    <- nrow(df_CHS_deterministic[[2]])
disp_lim    <- nrow(rc_bray_CHS %>% filter(RCbray > 0.95))
homo_disp   <- nrow(rc_bray_CHS %>% filter(RCbray < -0.95))
drift       <- nrow(rc_bray_CHS %>% filter(RCbray >= -0.95 & RCbray <= 0.95))

#repeat for each bioclimatic zone

# Build data with hierarchy + order
plot_data <- data.frame(
  
  Process = c(
    "Deterministic",
    "Stochastic",
    
    "Heterogeneous Selection",
    "Homogeneous Selection",
    
    "Dispersal Limitation",
    "Homogenizing Dispersal",
    "Drift"
  ),
  
  Size = c(
    det_hetero + det_homo,
    disp_lim + homo_disp + drift,
    
    det_hetero,
    det_homo,
    
    disp_lim,
    homo_disp,
    drift
  ),
  
  Ring = c(
    "Outer","Outer",
    "Inner","Inner",
    "Inner","Inner","Inner"
  ),
  
  Parent = c(
    NA, NA,
    "Deterministic","Deterministic",
    "Stochastic","Stochastic","Stochastic"
  ),
  
  Order = c(
    1, 2,   # Outer
    1, 2,   # Inside Deterministic
    3, 4, 5 # Inside Stochastic
  )
)

# Set factor levels explicitly
plot_data <- plot_data %>%
  arrange(Order) %>%
  mutate(
    Process = factor(Process, levels = Process)
  )

# Radius
plot_data$radius <- ifelse(plot_data$Ring == "Outer", 2, 1)

# Colors
cols <- c(
  "Deterministic"           = "#E74C3C",
  "Stochastic"              = "#3498DB",
  
  "Heterogeneous Selection" = "#E67E22",
  "Homogeneous Selection"   = "#F1C40F",
  
  "Dispersal Limitation"    = "#9B59B6",
  "Homogenizing Dispersal"  = "#2E86C1",
  "Drift"                   = "#1ABC9C"
)

# Plot
eco_process_CHS <- ggplot(plot_data,
                          aes(x = radius, y = Size, fill = Process)) +
  
  geom_col(width = 1, color = "white") +
  
  coord_polar(theta = "y") +
  
  xlim(0.5, 2.5) +
  
  scale_fill_manual(values = cols) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  
  ggtitle("Ecological processes for ecosystems / CHS")

eco_process_all/(eco_process_CCT+eco_process_CHS+eco_process_TEMP)




