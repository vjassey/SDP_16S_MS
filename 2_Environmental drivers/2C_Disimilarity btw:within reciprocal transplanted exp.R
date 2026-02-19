#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 2C - Disimilarity btw/within reciprocal transplanted exp
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Import data
load("RTExp_16s.RData")

###
### Prepare data frame for the plot
###
library(dplyr)
library(magrittr)
library(ggplot2)
library(microeco)
library(patchwork)


#hellinger transformation of the asv matrix
RTExp_16s_hel <- RTExp_16s
RTExp_16s_hel@otu_table <- otu_table(decostand(RTExp_16s_hel@otu_table, 'hellinger'), TRUE)

###
#### PCAO Based on relative abundances ----
###

# run pcoa
# make the PCoA using phyloseq
pcoa <- ordinate(RTExp_16s_hel, 'PCoA', 'horn')
pcoa_ordi12 <- plot_ordination(RTExp_16s_hel, pcoa, 'site')
pcoa_ordi13 <- plot_ordination(RTExp_16s_hel, pcoa, 'site', axes = c(4,5)) # run pcoa with Morisita-Horn distance

# get relative eigenvalues (% variance associated to different axes)
expl1 <- pcoa$values$Relative_eig[1] #axis1
expl2 <- pcoa$values$Relative_eig[2] #axis2
expl3 <- pcoa$values$Relative_eig[3] #axis3

#Extract the eigenvalues (variance explained) from the PCoA object
eigvals <- pcoa$values$Eigenvalues

# Calculate the proportion of variance explained by each axis
prop_var <- eigvals / sum(eigvals)

# Create a data frame for the first 20 axes
var_exp <- data.frame(
  Axis = 1:144,
  Eigenvalue = eigvals[1:144],
  Proportion = prop_var[1:144],
  Cumulative = cumsum(prop_var[1:144])
)

# Print the table
print(var_exp)

# Optionally, plot the variance explained
plot(var_exp$Axis, var_exp$Cumulative*100, type = "b", xlab = "Axis", ylab = "Cumulative Variance Explained",
     main = "Variance Explained by First 20 PCoA Axes")

## Target 100 axes to catch relevant centroids based on 80% variance explained! 

###
#### Distance between centroids from PCoA ----
###

# Combine PCoA results with treatment information
pcoa_data <- data.frame(
  pcoa$vectors[,1:100],
  Treatment = sample_data(RTExp_16s_hel)$Treatment,
  Destination = sample_data(RTExp_16s_hel)$Destination,
  Origin = sample_data(RTExp_16s_hel)$Origin)

# Calculate centroids for each treatment
#centroids <- aggregate(
#  . ~ Treatment+Destination+Origin,
#  data = pcoa_data[, c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5","PCoA6","PCoA7","PCoA8","PCoA9","PCoA10","PCoA11","PCoA12","Treatment", "Destination", "Origin")],
#  FUN = mean
#)

# Select the first 100 columns (assuming they are numeric)
numeric_cols <- pcoa_data[, 1:100]
# Select the grouping columns
grouping_cols <- pcoa_data[, c("Treatment", "Destination", "Origin")]
# Combine the numeric and grouping columns for aggregation
data_to_aggregate <- cbind(grouping_cols, numeric_cols)
# Aggregate using the mean function
centroids <- aggregate(
  . ~ Treatment + Destination + Origin,
  data = data_to_aggregate,
  FUN = mean
)


# Function to calculate Euclidean distance between two points
euclidean_distance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# Calculate the distance matrix using the dist() function
distance_matrix <- dist(as.matrix(centroids[, 4:103]), method = "euclidean")
# Convert the distance matrix to a square matrix
distance_matrix <- as.matrix(distance_matrix)
colnames(distance_matrix) <- centroids$Treatment
rownames(distance_matrix) <- centroids$Treatment

# Print the distance matrix
head(distance_matrix)

# Convert the distance matrix to a long-format data frame
pairwise_distances <- as.data.frame(as.table(distance_matrix))
names(pairwise_distances) <- c("Treatment1", "Treatment2", "Distance")


# Print the resulting data frame
head(pairwise_distances)


# Make comparison label like "A2C_C2C"
pairwise_distances <- pairwise_distances %>%
  mutate(comp = paste0(Treatment1, "_", Treatment2))

df <- pairwise_distances[,c(4,3)]

# Split Treatment_comp
df <- df %>% separate(comp, into = c("T1", "T2"), sep = "_")

###
#### Comparison Donor transplanted/Receptor versus Donor Transplanted/Donor Origin
###

# Manually define your pairings (comp1 and comp2) based on your example
pairings <- data.frame(
  comp1 = c("A2C_C2C","K2C_C2C","M2C_C2C","S2C_C2C",
            "A2K_K2K","C2K_K2K","M2K_K2K","S2K_K2K",
            "A2M_M2M", "C2M_M2M", "K2M_M2M", "S2M_M2M",
            "C2A_A2A", "K2A_A2A", "M2A_A2A", "S2A_A2A",
            "C2S_S2S", "K2S_S2S", "M2S_S2S", "A2S_S2S"),
  comp2 = c("A2C_A2A","K2C_K2K","M2C_M2M","S2C_S2S",
            "A2K_A2A","C2K_C2C","M2K_M2M","S2K_S2S",
            "A2M_A2A", "C2M_C2C", "K2M_K2K", "S2M_S2S",
            "C2A_C2C", "K2A_K2K", "M2A_M2M", "S2A_S2S",
            "C2S_C2C", "K2S_K2K", "M2S_M2M", "A2S_A2A"),
  stringsAsFactors = FALSE
)


# Join with df to get D1
pairings <- pairings %>%
  left_join(df %>% mutate(comp1_full = paste(T1,T2,sep="_")) %>% select(comp1_full, D1 = Distance),
            by = c("comp1" = "comp1_full")) %>%
  # Join with df to get D2
  left_join(df %>% mutate(comp2_full = paste(T1,T2,sep="_")) %>% select(comp2_full, D2 = Distance),
            by = c("comp2" = "comp2_full"))


head(pairings)

# Rename the distance columns to a common name
pairings_comp1 <- pairings[,c(1,3)] %>% dplyr::rename(Comp = comp1, Distance = D1)
pairings_comp1$Type <- 'Transplanted_vs_Receptor'
pairings_comp2 <- pairings[,c(2,4)] %>% dplyr::rename(Comp = comp2, Distance = D2)
pairings_comp2$Type <- 'Transplanted_vs_Origin'

# Combine the data frames
pairings_combined <- rbind(pairings_comp2, pairings_comp1)
str(pairings_combined)

boxplot(pairings_combined$Distance~as.factor(pairings_combined$Type))
summary(aov(pairings_combined$Distance~as.factor(pairings_combined$Type))) #pvalue < 0.004***
Anova(aov(pairings_combined$Distance~as.factor(pairings_combined$Type)))

# Statistical test
library(rstatix)
stat.test <- pairings_combined %>%
  anova_test(Distance~Type) %>%
  add_significance()
stat.test

#Plot
pcoa_comp1 <- ggplot(pairings_combined, aes(x = Type, y = Distance)) +
  geom_jitter(aes(col = Type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = Type)) +
  scale_color_manual(values = c("#2297E6", "#F5C710")) +
  scale_fill_manual(values = c("#2297E6", "#F5C710")) +
  theme_classic() +
  ylab(bquote(bold("Distance between centroids")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Receptor"))
pcoa_comp1  <- pcoa_comp1 + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                               label = "p.signif", 
                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                  symbols = c("****", "***", "**", "*", "ns")))


#IQR values
mean(pairings_combined$Distance[1:20])
quantile(pairings_combined$Distance[1:20], c(0.25))
quantile(pairings_combined$Distance[1:20], c(0.75))


###
#### Comparison with other transplanted communities from other donor sites
###

# Create a data frame for comp2 and join with df to get D2
pairings2_comp2 <- data.frame(
  comp = c("A2C_A2A","K2C_K2K","M2C_M2M","S2C_S2S",
           "A2K_A2A","C2K_C2C","M2K_M2M","S2K_S2S",
           "A2M_A2A", "C2M_C2C", "K2M_K2K", "S2M_S2S",
           "C2A_C2C", "K2A_K2K", "M2A_M2M", "S2A_S2S",
           "C2S_C2C", "K2S_K2K", "M2S_M2M", "A2S_A2A"),
  type = "Transplanted_vs_Origin",
  stringsAsFactors = FALSE
)

pairings2_comp2 <- pairings2_comp2 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D2 = Distance),
            by = c("comp" = "comp_full"))

# Create a data frame for comp3 and join with df to get D3
pairings2_comp3 <- data.frame(
  comp = c("A2C_K2C","A2C_M2C","A2C_S2C","A2C_C2C",
           "K2C_C2C","K2C_M2C","K2C_S2C","M2C_C2C",
           "M2C_S2C","S2C_C2C",
           "A2K_K2K","A2K_M2K","A2K_S2K","A2K_C2K",
           "C2K_K2K","C2K_M2K","C2K_S2K","M2K_K2K",
           "M2K_S2K","S2K_C2K",
           "A2M_M2M","A2M_C2M","A2M_K2M","A2M_S2M",
           "C2M_M2M","C2M_K2M","C2M_S2M","K2M_M2M",
           "K2M_S2M","S2M_M2M",
           "A2S_S2S","A2S_C2S","A2S_K2S","A2S_S2S",
           "C2S_S2S","C2S_K2S","C2S_S2S","K2S_M2S",
           "K2S_S2S","M2S_S2S",
           "C2A_A2A","C2A_K2A","C2A_M2A","C2A_S2A",
           "K2A_A2A","K2A_M2A","K2A_S2A","M2A_S2A",
           "M2A_A2A","S2A_A2A"),
  type = "Transplanted_vs_Transplanted",
  stringsAsFactors = FALSE
)

pairings2_comp3 <- pairings2_comp3 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D3 = Distance),
            by = c("comp" = "comp_full"))

# Rename the distance columns to a common name
pairings2_comp2 <- pairings2_comp2 %>% dplyr::rename(Distance = D2)
pairings2_comp3 <- pairings2_comp3 %>% dplyr::rename(Distance = D3)

# Combine the data frames
pairings2_combined <- rbind(pairings2_comp2, pairings2_comp3)

# View the result
head(pairings2_combined)

summary(aov(pairings2_combined$Distance~pairings2_combined$type))
Anova(aov(pairings2_combined$Distance~pairings2_combined$type), type="III")

# Statistical test
library(rstatix)
stat.test <- pairings2_combined %>%
  anova_test(Distance~type) %>%
  add_significance()
stat.test

# Extract the p-value for annotation
stars <- stat.test$p.signif

#Plot
pcoa_comp2 <- ggplot(pairings2_combined, aes(x = type, y = Distance)) +
  geom_jitter(aes(col = type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = type)) +
  scale_color_manual(values = c("#2297E6", "#DF536B")) +
  scale_fill_manual(values = c("#2297E6", "#DF536B")) +
  theme_classic() +
  ylab(bquote(bold("Distance between centroids")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Transplanted"))
pcoa_comp2  <- pcoa_comp2 + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                               label = "p.signif", 
                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                  symbols = c("****", "***", "**", "*", "ns")))



###
#### PCOA Based on P/A ----
###

# run pcoa
# make the PCoA using phyloseq
pcoa <- ordinate(RTExp_16s_hel, 'PCoA', 'jaccard')
pcoa_ordi12 <- plot_ordination(RTExp_16s_hel, pcoa, 'site')
pcoa_ordi13 <- plot_ordination(RTExp_16s_hel, pcoa, 'site', axes = c(4,5)) # run pcoa with Morisita-Horn distance

# get relative eigenvalues (% variance associated to different axes)
expl1 <- pcoa$values$Relative_eig[1] #axis1
expl2 <- pcoa$values$Relative_eig[2] #axis2
expl3 <- pcoa$values$Relative_eig[3] #axis3

#Extract the eigenvalues (variance explained) from the PCoA object
eigvals <- pcoa$values$Eigenvalues

# Calculate the proportion of variance explained by each axis
prop_var <- eigvals / sum(eigvals)

# Create a data frame for the first 20 axes
var_exp <- data.frame(
  Axis = 1:144,
  Eigenvalue = eigvals[1:144],
  Proportion = prop_var[1:144],
  Cumulative = cumsum(prop_var[1:144])
)

# Print the table
print(var_exp)

# Optionally, plot the variance explained
plot(var_exp$Axis, var_exp$Cumulative*100, type = "b", xlab = "Axis", ylab = "Cumulative Variance Explained",
     main = "Variance Explained by First 20 PCoA Axes")

## Target 100 axes to catch relevant centroids based on 80% variance explained! 

###
#### Distance between centroids from PCoA ----
###

# Combine PCoA results with treatment information
pcoa_data_pa <- data.frame(
  pcoa$vectors[,1:110],
  Treatment = sample_data(RTExp_16s_hel)$Treatment,
  Destination = sample_data(RTExp_16s_hel)$Destination,
  Origin = sample_data(RTExp_16s_hel)$Origin)

# Calculate centroids for each treatment
#centroids <- aggregate(
#  . ~ Treatment+Destination+Origin,
#  data = pcoa_data[, c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5","PCoA6","PCoA7","PCoA8","PCoA9","PCoA10","PCoA11","PCoA12","Treatment", "Destination", "Origin")],
#  FUN = mean
#)

# Select the first 100 columns (assuming they are numeric)
numeric_cols <- pcoa_data_pa[, 1:110]
# Select the grouping columns
grouping_cols <- pcoa_data_pa[, c("Treatment", "Destination", "Origin")]
# Combine the numeric and grouping columns for aggregation
data_to_aggregate <- cbind(grouping_cols, numeric_cols)
# Aggregate using the mean function
centroids_pa <- aggregate(
  . ~ Treatment + Destination + Origin,
  data = data_to_aggregate,
  FUN = mean
)
head(centroids_pa)

# Calculate the distance matrix using the dist() function
distance_matrix_pa <- dist(as.matrix(centroids_pa[, 4:113]), method = "euclidean")
# Convert the distance matrix to a square matrix
distance_matrix_pa <- as.matrix(distance_matrix_pa)
colnames(distance_matrix_pa) <- centroids_pa$Treatment
rownames(distance_matrix_pa) <- centroids_pa$Treatment

# Print the distance matrix
head(distance_matrix_pa)

# Convert the distance matrix to a long-format data frame
pairwise_distances_pa <- as.data.frame(as.table(distance_matrix_pa))
names(pairwise_distances_pa) <- c("Treatment1", "Treatment2", "Distance")


# Print the resulting data frame
head(pairwise_distances_pa)


# Make comparison label like "A2C_C2C"
pairwise_distances_pa <- pairwise_distances_pa %>%
  mutate(comp = paste0(Treatment1, "_", Treatment2))

df <- pairwise_distances_pa[,c(4,3)]

# Split Treatment_comp
df <- df %>% separate(comp, into = c("T1", "T2"), sep = "_")

###
#### Comparison Donor transplanted/Receptor versus Donor Transplanted/Donor Origin
###

# Manually define your pairings (comp1 and comp2) based on your example
pairings <- data.frame(
  comp1 = c("A2C_C2C","K2C_C2C","M2C_C2C","S2C_C2C",
            "A2K_K2K","C2K_K2K","M2K_K2K","S2K_K2K",
            "A2M_M2M", "C2M_M2M", "K2M_M2M", "S2M_M2M",
            "C2A_A2A", "K2A_A2A", "M2A_A2A", "S2A_A2A",
            "C2S_S2S", "K2S_S2S", "M2S_S2S", "A2S_S2S"),
  comp2 = c("A2C_A2A","K2C_K2K","M2C_M2M","S2C_S2S",
            "A2K_A2A","C2K_C2C","M2K_M2M","S2K_S2S",
            "A2M_A2A", "C2M_C2C", "K2M_K2K", "S2M_S2S",
            "C2A_C2C", "K2A_K2K", "M2A_M2M", "S2A_S2S",
            "C2S_C2C", "K2S_K2K", "M2S_M2M", "A2S_A2A"),
  stringsAsFactors = FALSE
)


# Join with df to get D1
pairings <- pairings %>%
  left_join(df %>% mutate(comp1_full = paste(T1,T2,sep="_")) %>% select(comp1_full, D1 = Distance),
            by = c("comp1" = "comp1_full")) %>%
  # Join with df to get D2
  left_join(df %>% mutate(comp2_full = paste(T1,T2,sep="_")) %>% select(comp2_full, D2 = Distance),
            by = c("comp2" = "comp2_full"))


head(pairings)

# Rename the distance columns to a common name
pairings_comp1 <- pairings[,c(1,3)] %>% dplyr::rename(Comp = comp1, Distance = D1)
pairings_comp1$Type <- 'Transplanted_vs_Receptor'
pairings_comp2 <- pairings[,c(2,4)] %>% dplyr::rename(Comp = comp2, Distance = D2)
pairings_comp2$Type <- 'Transplanted_vs_Origin'

# Combine the data frames
pairings_combined <- rbind(pairings_comp2, pairings_comp1)
str(pairings_combined)

boxplot(pairings_combined$Distance~as.factor(pairings_combined$Type))
summary(aov(pairings_combined$Distance~as.factor(pairings_combined$Type))) #pvalue < 0.09
Anova(aov(pairings_combined$Distance~as.factor(pairings_combined$Type)))

# Statistical test
library(rstatix)
stat.test <- pairings_combined %>%
  anova_test(Distance~Type) %>%
  add_significance()
stat.test

#Plot
pcoa_comp1_pa <- ggplot(pairings_combined, aes(x = Type, y = Distance)) +
  geom_jitter(aes(col = Type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = Type)) +
  scale_color_manual(values = c("#2297E6", "#F5C710")) +
  scale_fill_manual(values = c("#2297E6", "#F5C710")) +
  theme_classic() +
  ylab(bquote(bold("Distance between centroids")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Receptor"))
pcoa_comp1_pa  <- pcoa_comp1_pa + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                                     label = "p.signif", 
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))

###
#### Comparison with other transplanted communities from other donor sites
###

# Create a data frame for comp2 and join with df to get D2
pairings2_comp2 <- data.frame(
  comp = c("A2C_A2A","K2C_K2K","M2C_M2M","S2C_S2S",
           "A2K_A2A","C2K_C2C","M2K_M2M","S2K_S2S",
           "A2M_A2A", "C2M_C2C", "K2M_K2K", "S2M_S2S",
           "C2A_C2C", "K2A_K2K", "M2A_M2M", "S2A_S2S",
           "C2S_C2C", "K2S_K2K", "M2S_M2M", "A2S_A2A"),
  type = "Transplanted_vs_Origin",
  stringsAsFactors = FALSE
)

pairings2_comp2 <- pairings2_comp2 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D2 = Distance),
            by = c("comp" = "comp_full"))

# Create a data frame for comp3 and join with df to get D3
pairings2_comp3 <- data.frame(
  comp = c("A2C_K2C","A2C_M2C","A2C_S2C","A2C_C2C",
           "K2C_C2C","K2C_M2C","K2C_S2C","M2C_C2C",
           "M2C_S2C","S2C_C2C",
           "A2K_K2K","A2K_M2K","A2K_S2K","A2K_C2K",
           "C2K_K2K","C2K_M2K","C2K_S2K","M2K_K2K",
           "M2K_S2K","S2K_C2K",
           "A2M_M2M","A2M_C2M","A2M_K2M","A2M_S2M",
           "C2M_M2M","C2M_K2M","C2M_S2M","K2M_M2M",
           "K2M_S2M","S2M_M2M",
           "A2S_S2S","A2S_C2S","A2S_K2S","A2S_S2S",
           "C2S_S2S","C2S_K2S","C2S_S2S","K2S_M2S",
           "K2S_S2S","M2S_S2S",
           "C2A_A2A","C2A_K2A","C2A_M2A","C2A_S2A",
           "K2A_A2A","K2A_M2A","K2A_S2A","M2A_S2A",
           "M2A_A2A","S2A_A2A"),
  type = "Transplanted_vs_Transplanted",
  stringsAsFactors = FALSE
)

pairings2_comp3 <- pairings2_comp3 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D3 = Distance),
            by = c("comp" = "comp_full"))

# Rename the distance columns to a common name
pairings2_comp2 <- pairings2_comp2 %>% dplyr::rename(Distance = D2)
pairings2_comp3 <- pairings2_comp3 %>% dplyr::rename(Distance = D3)

# Combine the data frames
pairings2_combined <- rbind(pairings2_comp2, pairings2_comp3)

# View the result
head(pairings2_combined)

summary(aov(pairings2_combined$Distance~pairings2_combined$type))
Anova(aov(pairings2_combined$Distance~pairings2_combined$type), type="III")

# Statistical test
library(rstatix)
stat.test <- pairings2_combined %>%
  anova_test(Distance~type) %>%
  add_significance()
stat.test

# Extract the p-value for annotation
stars <- stat.test$p.signif

#Plot
pcoa_comp2_pa <- ggplot(pairings2_combined, aes(x = type, y = Distance)) +
  geom_jitter(aes(col = type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = type)) +
  scale_color_manual(values = c("#2297E6", "#DF536B")) +
  scale_fill_manual(values = c("#2297E6", "#DF536B")) +
  theme_classic() +
  ylab(bquote(bold("Distance between centroids")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Transplanted"))
pcoa_comp2_pa  <- pcoa_comp2_pa + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                                     label = "p.signif", 
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))


#build plot
pcoa_comp1+pcoa_comp2+pcoa_comp1_pa+pcoa_comp2_pa


###
#### Pairwise Betadiversity Comparisons ----
###

# create df and merge with sample info
df_ordihorn <- as.data.frame(as.matrix(sample_data(RTExp_16s_hel)))
df_ordihorn$Treatment <- as.factor(df_ordihorn$Treatment)

# Calculate pairwise beta-diversity (e.g., Bray-Curtis)
beta_diversity <- vegdist(t(RTExp_16s_hel@otu_table), method = "horn")

# Convert the dissimilarity matrix to a data frame
beta_df <- as.matrix(beta_diversity)


# Create a data frame for plotting
treatments <- df_ordihorn$Treatment
plot_data <- data.frame(
  Sample1 = rep(rownames(beta_df), each = ncol(beta_df)),
  Sample2 = rep(colnames(beta_df), times = nrow(beta_df)),
  Dissimilarity = as.vector(beta_df),
  Treatment1 = rep(treatments, each = ncol(beta_df)),
  Treatment2 = rep(treatments, times = nrow(beta_df))
)
#write.table(plot_data, "beta-div_mixopeat.csv", sep = ";")


# View the result
#head(aggregated_dissimilarity)

df <- plot_data[,c(3,1,2)]
colnames(df) <- c("Distance", "T1", "T2")

###
#### Comparison Donor transplanted/Receptor versus Donor Transplanted/Donor Origin
###

# Manually define your pairings (comp1 and comp2) based on your example
pairings <- data.frame(
  comp1 = c("A2C1_C2C1","K2C1_C2C1","M2C1_C2C1","S2C1_C2C1",
            "A2K1_K2K1","C2K1_K2K1","M2K1_K2K1","S2K1_K2K1",
            "A2M1_M2M1", "C2M1_M2M1", "K2M1_M2M1", "S2M1_M2M1",
            "C2A1_A2A1", "K2A1_A2A1", "M2A1_A2A1", "S2A1_A2A1",
            "C2S1_S2S1", "K2S1_S2S1", "M2S1_S2S1", "A2S1_S2S1",
            "A2C2_C2C2","K2C2_C2C2","M2C2_C2C2","S2C2_C2C2",
            "A2K2_K2K2","C2K2_K2K2","M2K2_K2K2","S2K2_K2K2",
            "A2M2_M2M2", "C2M2_M2M2", "K2M2_M2M2", "S2M2_M2M2",
            "C2A2_A2A2", "K2A2_A2A2", "M2A2_A2A2", "S2A2_A2A2",
            "C2S2_S2S2", "K2S2_S2S2", "M2S2_S2S2", "A2S2_S2S2",
            "A2C3_C2C3","K2C3_C2C3","M2C3_C2C3","S2C3_C2C3",
            "A2K3_K2K3","C2K3_K2K3","M2K3_K2K3","S2K3_K2K3",
            "A2M3_M2M3", "C2M3_M2M3", "K2M3_M2M3", "S2M3_M2M3",
            "C2A3_A2A3", "K2A3_A2A3", "M2A3_A2A3", "S2A3_A2A3",
            "C2S3_S2S3", "K2S3_S2S3", "M2S3_S2S3", "A2S3_S2S3",
            "A2C4_C2C4","K2C4_C2C4","M2C4_C2C4","S2C4_C2C4",
            "A2K4_K2K4","C2K4_K2K4","M2K4_K2K4","S2K4_K2K4",
            "A2M4_M2M4", "C2M4_M2M4", "K2M4_M2M4", "S2M4_M2M4",
            "C2A4_A2A4", "K2A4_A2A4", "M2A4_A2A4", "S2A4_A2A4",
            "C2S4_S2S4", "K2S4_S2S4", "M2S4_S2S4", "A2S4_S2S4",
            "A2C5_C2C5","K2C5_C2C5","M2C5_C2C5","S2C5_C2C5",
            "A2K5_K2K5","C2K5_K2K5","M2K5_K2K5","S2K5_K2K5",
            "A2M5_M2M5", "C2M5_M2M5", "K2M5_M2M5", "S2M5_M2M5",
            "C2A5_A2A5", "K2A5_A2A5", "M2A5_A2A5", "S2A5_A2A5",
            "C2S5_S2S5", "K2S5_S2S5", "M2S5_S2S5", "A2S5_S2S5"),
  comp2 = c("A2C1_A2A1","K2C1_K2K1","M2C1_M2M1","S2C1_S2S1",
            "A2K1_A2A1","C2K1_C2C1","M2K1_M2M1","S2K1_S2S1",
            "A2M1_A2A1", "C2M1_C2C1", "K2M1_K2K1", "S2M1_S2S1",
            "C2A1_C2C1", "K2A1_K2K1", "M2A1_M2M1", "S2A1_S2S1",
            "C2S2_C2C2", "K2S2_K2K2", "M2S2_M2M2", "A2S2_A2A2",
            "A2C2_A2A2","K2C2_K2K2","M2C2_M2M2","S2C2_S2S2",
            "A2K2_A2A2","C2K2_C2C2","M2K2_M2M2","S2K2_S2S2",
            "A2M2_A2A2", "C2M2_C2C2", "K2M2_K2K2", "S2M2_S2S2",
            "C2A2_C2C2", "K2A2_K2K2", "M2A2_M2M2", "S2A2_S2S2",
            "C2S2_C2C2", "K2S2_K2K2", "M2S2_M2M2", "A2S2_A2A2",
            "A2C3_A2A3","K2C3_K2K3","M2C3_M2M3","S2C3_S2S3",
            "A2K3_A2A3","C2K3_C2C3","M2K3_M2M3","S2K3_S2S3",
            "A2M3_A2A3", "C2M3_C2C3", "K2M3_K2K3", "S2M3_S2S3",
            "C2A3_C2C3", "K2A3_K2K3", "M2A3_M2M3", "S2A3_S2S3",
            "C2S3_C2C3", "K2S3_K2K3", "M2S3_M2M3", "A2S3_A2A3",
            "A2C4_A2A4","K2C4_K2K4","M2C4_M2M4","S2C4_S2S4",
            "A2K4_A2A4","C2K4_C2C4","M2K4_M2M4","S2K4_S2S4",
            "A2M4_A2A4", "C2M4_C2C4", "K2M4_K2K4", "S2M4_S2S4",
            "C2A4_C2C4", "K2A4_K2K4", "M2A4_M2M4", "S2A4_S2S4",
            "C2S4_C2C4", "K2S4_K2K4", "M2S4_M2M4", "A2S4_A2A4",
            "A2C5_A2A5","K2C5_K2K5","M2C5_M2M5","S2C5_S2S5",
            "A2K5_A2A5","C2K5_C2C5","M2K5_M2M5","S2K5_S2S5",
            "A2M5_A2A5", "C2M5_C2C5", "K2M5_K2K5", "S2M5_S2S5",
            "C2A5_C2C5", "K2A5_K2K5", "M2A5_M2M5", "S2A5_S2S5",
            "C2S5_C2C5", "K2S5_K2K5", "M2S5_M2M5", "A2S5_A2A5"),
  stringsAsFactors = FALSE
)

head(pairings)

pairings <- pairings %>%
  left_join(df %>% mutate(comp1_full = paste(T1,T2,sep="_")) %>% select(comp1_full, D1 = Distance),
            by = c("comp1" = "comp1_full")) %>%
  # Join with df to get D2
  left_join(df %>% mutate(comp2_full = paste(T1,T2,sep="_")) %>% select(comp2_full, D2 = Distance),
            by = c("comp2" = "comp2_full"))

# Rename the distance columns to a common name
pairings_comp1 <- pairings[,c(1,3)] %>% dplyr::rename(Comp = comp1, Distance = D1)
pairings_comp1$Type <- 'Transplanted_vs_Receptor'
pairings_comp2 <- pairings[,c(2,4)] %>% dplyr::rename(Comp = comp2, Distance = D2)
pairings_comp2$Type <- 'Transplanted_vs_Origin'

# Combine the data frames
pairings_combined <- rbind(pairings_comp2, pairings_comp1)
str(pairings_combined)

summary(aov(pairings_combined$Distance~as.factor(pairings_combined$Type))) #pvalue < 0.06

# Statistical test
library(rstatix)
stat.test <- pairings_combined %>%
  anova_test(Distance~Type) %>%
  add_significance()
stat.test

#Plot
diss_comp1 <- ggplot(pairings_combined, aes(x = Type, y = Distance)) +
  geom_jitter(aes(col = Type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = Type)) +
  scale_color_manual(values = c("#2297E6", "#F5C710")) +
  scale_fill_manual(values = c("#2297E6", "#F5C710")) +
  theme_classic() +
  ylab(bquote(bold("Jaccard's distance")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Receptor"))
diss_comp1  <- diss_comp1 + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                               label = "p.signif", 
                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                  symbols = c("****", "***", "**", "*", "ns")))
diss_comp1

diss_comp1_pa <- ggplot(pairings_combined, aes(x = Type, y = Distance)) +
  geom_jitter(aes(col = Type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = Type)) +
  scale_color_manual(values = c("#2297E6", "#F5C710")) +
  scale_fill_manual(values = c("#2297E6", "#F5C710")) +
  theme_classic() +
  ylab(bquote(bold("Jaccard's dissimilarity")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Receptor"))
diss_comp1_pa  <- diss_comp1_pa + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                                     label = "p.signif", 
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))
diss_comp1_pa 

###
#### Comparison with other transplanted communities from other donor sites
###

# Create a data frame for comp2 and join with df to get D2
pairings2_comp2 <- data.frame(
  comp = c("A2C1_A2A1","K2C1_K2K1","M2C1_M2M1","S2C1_S2S1",
           "A2K1_A2A1","C2K1_C2C1","M2K1_M2M1","S2K1_S2S1",
           "A2M1_A2A1", "C2M1_C2C1", "K2M1_K2K1", "S2M1_S2S1",
           "C2A1_C2C1", "K2A1_K2K1", "M2A1_M2M1", "S2A1_S2S1",
           "C2S2_C2C2", "K2S2_K2K2", "M2S2_M2M2", "A2S2_A2A2",
           "A2C2_A2A2","K2C2_K2K2","M2C2_M2M2","S2C2_S2S2",
           "A2K2_A2A2","C2K2_C2C2","M2K2_M2M2","S2K2_S2S2",
           "A2M2_A2A2", "C2M2_C2C2", "K2M2_K2K2", "S2M2_S2S2",
           "C2A2_C2C2", "K2A2_K2K2", "M2A2_M2M2", "S2A2_S2S2",
           "C2S2_C2C2", "K2S2_K2K2", "M2S2_M2M2", "A2S2_A2A2",
           "A2C3_A2A3","K2C3_K2K3","M2C3_M2M3","S2C3_S2S3",
           "A2K3_A2A3","C2K3_C2C3","M2K3_M2M3","S2K3_S2S3",
           "A2M3_A2A3", "C2M3_C2C3", "K2M3_K2K3", "S2M3_S2S3",
           "C2A3_C2C3", "K2A3_K2K3", "M2A3_M2M3", "S2A3_S2S3",
           "C2S3_C2C3", "K2S3_K2K3", "M2S3_M2M3", "A2S3_A2A3",
           "A2C4_A2A4","K2C4_K2K4","M2C4_M2M4","S2C4_S2S4",
           "A2K4_A2A4","C2K4_C2C4","M2K4_M2M4","S2K4_S2S4",
           "A2M4_A2A4", "C2M4_C2C4", "K2M4_K2K4", "S2M4_S2S4",
           "C2A4_C2C4", "K2A4_K2K4", "M2A4_M2M4", "S2A4_S2S4",
           "C2S4_C2C4", "K2S4_K2K4", "M2S4_M2M4", "A2S4_A2A4",
           "A2C5_A2A5","K2C5_K2K5","M2C5_M2M5","S2C5_S2S5",
           "A2K5_A2A5","C2K5_C2C5","M2K5_M2M5","S2K5_S2S5",
           "A2M5_A2A5", "C2M5_C2C5", "K2M5_K2K5", "S2M5_S2S5",
           "C2A5_C2C5", "K2A5_K2K5", "M2A5_M2M5", "S2A5_S2S5",
           "C2S5_C2C5", "K2S5_K2K5", "M2S5_M2M5", "A2S5_A2A5"),
  type = "Transplanted_vs_Origin",
  stringsAsFactors = FALSE
)

pairings2_comp2 <- pairings2_comp2 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D2 = Distance),
            by = c("comp" = "comp_full"))

# Create a data frame for comp3 and join with df to get D3
pairings2_comp3 <- data.frame(
  comp = c("A2C1_K2C1","A2C1_M2C1","A2C1_S2C1","A2C1_C2C1",
           "K2C1_C2C1","K2C1_M2C1","K2C1_S2C1","M2C1_C2C1",
           "M2C1_S2C1","S2C1_C2C1",
           "A2K1_K2K1","A2K1_M2K1","A2K1_S2K1","A2K1_C2K1",
           "C2K1_K2K1","C2K1_M2K1","C2K1_S2K1","M2K1_K2K1",
           "M2K1_S2K1","S2K1_C2K1",
           "A2M1_M2M1","A2M1_C2M1","A2M1_K2M1","A2M1_S2M1",
           "C2M1_M2M1","C2M1_K2M1","C2M1_S2M1","K2M1_M2M1",
           "K2M1_S2M1","S2M1_M2M1",
           "A2S1_S2S1","A2S1_C2S1","A2S1_K2S1","A2S1_S2S1",
           "C2S1_S2S1","C2S1_K2S1","C2S1_S2S1","K2S1_M2S1",
           "K2S1_S2S1","M2S1_S2S1",
           "C2A1_A2A1","C2A1_K2A1","C2A1_M2A1","C2A1_S2A1",
           "K2A1_A2A1","K2A1_M2A1","K2A1_S2A1","M2A1_S2A1",
           "M2A1_A2A1","S2A1_A2A1",
           "A2C2_K2C2","A2C2_M2C2","A2C2_S2C2","A2C2_C2C2",
           "K2C2_C2C2","K2C2_M2C2","K2C2_S2C2","M2C2_C2C2",
           "M2C2_S2C2","S2C2_C2C2",
           "A2K2_K2K2","A2K2_M2K2","A2K2_S2K2","A2K2_C2K2",
           "C2K2_K2K2","C2K2_M2K2","C2K2_S2K2","M2K2_K2K2",
           "M2K2_S2K2","S2K2_C2K2",
           "A2M2_M2M2","A2M2_C2M2","A2M2_K2M2","A2M2_S2M2",
           "C2M2_M2M2","C2M2_K2M2","C2M2_S2M2","K2M2_M2M2",
           "K2M2_S2M2","S2M2_M2M2",
           "A2S2_S2S2","A2S2_C2S2","A2S2_K2S2","A2S2_S2S2",
           "C2S2_S2S2","C2S2_K2S2","C2S2_S2S2","K2S2_M2S2",
           "K2S2_S2S2","M2S2_S2S2",
           "C2A2_A2A2","C2A2_K2A2","C2A2_M2A2","C2A2_S2A2",
           "K2A2_A2A2","K2A2_M2A2","K2A2_S2A2","M2A2_S2A2",
           "M2A2_A2A2","S2A2_A2A2",
           "A2C3_K2C3","A2C3_M2C3","A2C3_S2C3","A2C3_C2C3",
           "K2C3_C2C3","K2C3_M2C3","K2C3_S2C3","M2C3_C2C3",
           "M2C3_S2C3","S2C3_C2C3",
           "A2K3_K2K3","A2K3_M2K3","A2K3_S2K3","A2K3_C2K3",
           "C2K3_K2K3","C2K3_M2K3","C2K3_S2K3","M2K3_K2K3",
           "M2K3_S2K3","S2K3_C2K3",
           "A2M3_M2M3","A2M3_C2M3","A2M3_K2M3","A2M3_S2M3",
           "C2M3_M2M3","C2M3_K2M3","C2M3_S2M3","K2M3_M2M3",
           "K2M3_S2M3","S2M3_M2M3",
           "A2S3_S2S3","A2S3_C2S3","A2S3_K2S3","A2S3_S2S3",
           "C2S3_S2S3","C2S3_K2S3","C2S3_S2S3","K2S3_M2S3",
           "K2S3_S2S3","M2S3_S2S3",
           "C2A3_A2A3","C2A3_K2A3","C2A3_M2A3","C2A3_S2A3",
           "K2A3_A2A3","K2A3_M2A3","K2A3_S2A3","M2A3_S2A3",
           "M2A3_A2A3","S2A3_A2A3",
           "A2C4_K2C4","A2C4_M2C4","A2C4_S2C4","A2C4_C2C4",
           "K2C4_C2C4","K2C4_M2C4","K2C4_S2C4","M2C4_C2C4",
           "M2C4_S2C4","S2C4_C2C4",
           "A2K4_K2K4","A2K4_M2K4","A2K4_S2K4","A2K4_C2K4",
           "C2K4_K2K4","C2K4_M2K4","C2K4_S2K4","M2K4_K2K4",
           "M2K4_S2K4","S2K4_C2K4",
           "A2M4_M2M4","A2M4_C2M4","A2M4_K2M4","A2M4_S2M4",
           "C2M4_M2M4","C2M4_K2M4","C2M4_S2M4","K2M4_M2M4",
           "K2M4_S2M4","S2M4_M2M4",
           "A2S4_S2S4","A2S4_C2S4","A2S4_K2S4","A2S4_S2S4",
           "C2S4_S2S4","C2S4_K2S4","C2S4_S2S4","K2S4_M2S4",
           "K2S4_S2S4","M2S4_S2S4",
           "C2A4_A2A4","C2A4_K2A4","C2A4_M2A4","C2A4_S2A4",
           "K2A4_A2A4","K2A4_M2A4","K2A4_S2A4","M2A4_S2A4",
           "M2A4_A2A4","S2A4_A2A4",
           "A2C5_K2C5","A2C5_M2C5","A2C5_S2C5","A2C5_C2C5",
           "K2C5_C2C5","K2C5_M2C5","K2C5_S2C5","M2C5_C2C5",
           "M2C5_S2C5","S2C5_C2C5",
           "A2K5_K2K5","A2K5_M2K5","A2K5_S2K5","A2K5_C2K5",
           "C2K5_K2K5","C2K5_M2K5","C2K5_S2K5","M2K5_K2K5",
           "M2K5_S2K5","S2K5_C2K5",
           "A2M5_M2M5","A2M5_C2M5","A2M5_K2M5","A2M5_S2M5",
           "C2M5_M2M5","C2M5_K2M5","C2M5_S2M5","K2M5_M2M5",
           "K2M5_S2M5","S2M5_M2M5",
           "A2S5_S2S5","A2S5_C2S5","A2S5_K2S5","A2S5_S2S5",
           "C2S5_S2S5","C2S5_K2S5","C2S5_S2S5","K2S5_M2S5",
           "K2S5_S2S5","M2S5_S2S5",
           "C2A5_A2A5","C2A5_K2A5","C2A5_M2A5","C2A5_S2A5",
           "K2A5_A2A5","K2A5_M2A5","K2A5_S2A5","M2A5_S2A5",
           "M2A5_A2A5","S2A5_A2A5"),
  type = "Trnasplanted_vs_Transplanted",
  stringsAsFactors = FALSE
)

pairings2_comp3 <- pairings2_comp3 %>%
  left_join(df %>% mutate(comp_full = paste(T1, T2, sep = "_")) %>% select(comp_full, D3 = Distance),
            by = c("comp" = "comp_full"))

# Rename the distance columns to a common name
pairings2_comp2 <- pairings2_comp2 %>% dplyr::rename(Distance = D2)
pairings2_comp3 <- pairings2_comp3 %>% dplyr::rename(Distance = D3)

# Combine the data frames
pairings2_combined <- rbind(pairings2_comp2, pairings2_comp3)

# View the result
head(pairings2_combined)

boxplot(pairings2_combined$Distance~pairings2_combined$type)
summary(aov(pairings2_combined$Distance~pairings2_combined$type))
Anova(aov(pairings2_combined$Distance~pairings2_combined$type))

# Statistical test
library(rstatix)
stat.test <- pairings2_combined %>%
  anova_test(Distance~type) %>%
  add_significance()
stat.test

# Extract the p-value for annotation
stars <- stat.test$p.signif

#Plot
diss_comp2 <- ggplot(pairings2_combined, aes(x = type, y = Distance)) +
  geom_jitter(aes(col = type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = type)) +
  scale_color_manual(values = c("#2297E6", "#DF536B")) +
  scale_fill_manual(values = c("#2297E6", "#DF536B")) +
  theme_classic() +
  ylab(bquote(bold("Jaccard's dissimilarity")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Transplanted"))
diss_comp2  <- diss_comp2_pa + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                                  label = "p.signif", 
                                                  symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                     symbols = c("****", "***", "**", "*", "ns")))
diss_comp2

diss_comp2_pa <- ggplot(pairings2_combined, aes(x = type, y = Distance)) +
  geom_jitter(aes(col = type), width = .1, height = 0, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, aes(fill = type)) +
  scale_color_manual(values = c("#2297E6", "#DF536B")) +
  scale_fill_manual(values = c("#2297E6", "#DF536B")) +
  theme_classic() +
  ylab(bquote(bold("Jaccard's dissimilarity")))

cmpr <- list(c("Transplanted_vs_Origin","Transplanted_vs_Transplanted"))
diss_comp2_pa  <- diss_comp2_pa + stat_compare_means(comparisons = cmpr, tip.length=0.01,
                                                     label = "p.signif", 
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))
diss_comp2_pa

#build plot
diss_comp1+diss_comp2+diss_comp1_pa+diss_comp2_pa

