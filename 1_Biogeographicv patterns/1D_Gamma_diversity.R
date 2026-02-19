#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 1D - Gamma diversity
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## This code is adapted from Ezzat, L., Peter, H., Bourquin, M. et al. 
#Diversity and biogeography of the bacterial microbiome in glacier-fed streams. 
#Nature 637, 622–630 (2025). https://doi.org/10.1038/s41586-024-08313-z

library(phyloseq)
library(phyloseqCompanion)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(DescTools)
library(devtools)
library(iNEXT)

## Load 16S processed data
load("SDPs_16s.RData") # input file

merged_NOMIS_DEBLUR <- SDPs_16s

sample_lists <- list()
asv_lists <- list()
vec_lists <- list()

# Loop through each site and perform the subset and merge operations
for (site in sites) {
  subset_result <- subset_samples(merged_NOMIS_DEBLUR, koppen_clim == site)
  merged_result <- merge_samples(subset_result,"Sample")
  
  # ## create asv table
  asv_table_site <- otu_table(merged_result, taxa_are_rows=T)
  asv_table_site_t <- t(asv_table_site)
  asv_lists[[site]] <- asv_table_site_t
  #
  # ## calculate rowSums and
  asv_table_df <- asv_table_site_t
  sumrow_site <- unname(rowSums(asv_table_df>0))
  sort_site<- sort(sumrow_site, decreasing=T)
  vec_site <- sort_site[sort_site >0]
  vec_lists[[site]] <- vec_site
}

list_exped_all <- list(temperate=c(ncol(asv_lists$Temperate),vec_lists$Temperate),cct=c(ncol(asv_lists$`Continental cold/Tundra`),vec_lists$`Continental cold/Tundra`),
                       chs=c(ncol(asv_lists$`Continental hot summer`),vec_lists$`Continental hot summer`))


out_all_exped <- iNEXT(list_exped_all, q=0, datatype="incidence_freq", se=T, conf=0.95, nboot=99)

df <- fortify(out_all_exped, type =1)

df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method, 
                         c("Rarefaction", "Extrapolation"),
)

#df.asympote <- data.frame(y = c(24,3), Asymptote = c("temperate","cct","chs"))

ggplot(df, aes(x=x, y=y, colour=Assemblage)) + 
  #geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype= Method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=Assemblage, colour=NULL), alpha=0.2) +
  labs(x="Number of peatland site", y="Species diversity") +
  scale_fill_manual(values=c("#6EB689",  "#7C8AB5","#C45858")
  )+
  scale_color_manual(values=c("#6EB689",  "#7C8AB5","#C45858")
  )+
  scale_linetype_discrete(name ="Method")+
  theme_bw() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

inext_freq_results <- out_all_exped$AsyEst  
inext_freq_results$prop <- inext_freq_results$Observed/inext_freq_results$Estimator
inext_freq_results<-inext_freq_results[inext_freq_results$Diversity == 'Species richness',]
median_GD_freq <- inext_freq_results %>% 
  summarise(med = median(prop), 
            lower_quartile = quantile(prop, 0.25),
            median = quantile(prop, 0.5),
            upper_quartile = quantile(prop, 0.75))


