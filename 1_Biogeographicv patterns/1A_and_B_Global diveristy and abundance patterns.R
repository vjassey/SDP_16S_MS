#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 1A and B - Global diversity and abundance patterns of the SDPs
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(ggplot2)
library(dplyr)
library(patchwork)

#import data
bac_biom <- read.csv("bacterial_biomass_data.csv", sep = ',') #input file
bac_richn <- read.csv("bacterial_richness_data.csv", sep = ',') #input file

###
#### Density plot abundance and ricnhess ----
###
# Density plot in ggplot2
theme_set(new = theme_classic())
dens_ab <- ggplot(bac_biom, aes(bac_biomass)) +
  geom_density(color="gray30", fill="gray30") 
dens_rich <- ggplot(bac_richn, aes(n0)) +
  geom_density(color="gray30", fill="gray30") 

# Set climate zones color
climate_palette <- paletteer::paletteer_d("MetBrewer::Austria")[1:3]
names(climate_palette) <- c("Continental hot summer","Continental cold/Tundra","Temperate")

# Richness
p_rich <- bac_richn %>%
  ggplot(aes(x = koppen_clim, y = n0, fill = koppen_clim)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(color = "black", 
              size = 1.5, 
              alpha = 0.7, 
              width = 0.2) +
  coord_flip()+
  scale_fill_manual(values = climate_palette,name="Climatic zone",guide="none")
# Add marginal density histogram at top
p_den_rich <- ggplot(bac_richn, aes(x = n0)) +
  geom_density(fill = "gray30", alpha = 0.8, color = "black") +
  theme_classic(base_size = 14) +
  labs(y = "Frequency") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0)
  )
# Display plot
p_final_rich <- p_den_rich / p_rich + plot_layout(heights = c(1, 4))


# Abundance
p_ab <- bac_biom %>%
  ggplot(aes(x = koppen_clim, y = bac_biomass, fill = koppen_clim)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(color = "black", 
              size = 1.5, 
              alpha = 0.7, 
              width = 0.2) +
  coord_flip()+
  scale_fill_manual(values = climate_palette,name="Climatic zone",guide="none")
# Add marginal density histogram at top
p_den_ab <- ggplot(bac_biom, aes(x = bac_biomass)) +
  geom_density(fill = "gray30", alpha = 0.8, color = "black") +
  theme_classic(base_size = 14) +
  labs(y = "Frequency") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0)
  )
# Display plot
p_final_ab <- p_den_ab / p_ab + plot_layout(heights = c(1, 4))
p_final_ab

p_final_rich | p_final_ab
