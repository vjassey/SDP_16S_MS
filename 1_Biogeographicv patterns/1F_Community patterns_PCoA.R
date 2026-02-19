#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 1F - PCoA analyses
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load packages =================================================
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(pairwiseAdonis)

# Custom fct =======================================
StatCentSeg <- ggplot2::ggproto("StatCentSeg", Stat,
                                compute_group = function(data, scales, params,
                                                         cfun=median) {
                                  data$xend <- cfun(data$x)
                                  data$yend <- cfun(data$y)
                                  return(data)
                                },
                                required_aes = c("x", "y")
)
stat_centseg <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA, 
                         inherit.aes = TRUE, cfun=median, ...) {
  layer(
    stat = StatCentSeg, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, cfun = cfun, ...)
  )
}


agg.table.taxo <- function(tab, tax.lvl="genus", tax.table) {
  tax.table <- tax.table[match(rownames(tab),
                               rownames(tax.table)),]
  message(paste('Table aggregation to the', tax.lvl, "level."))
  #message('Please be sure that the ASV/OTU table and the taxonomy table are ordered the same way')
  if(nrow(tab) != nrow(tax.table)) stop("The ASV/OTU table and the taxonomy table do not have the same number of rows")
  tax <- tax.table[,grep(tax.lvl, colnames(tax.table), ignore.case = T)]
  tax[is.na(tax)] <- "Unknown"
  tab <- aggregate(tab, by=list("taxo"=tax), FUN=sum)
  rownames(tab) <- tab[,1]  
  tab <- tab[,-1]
  return(tab)
}


## Load 16S processed data =========================
load("SDPs_16s.RData") # input file


# Prepare for plot================================================
tax.lvl <- "Class_conf_rdp"

tab <- SDPs_16s@otu_table
class_tot <- agg.table.taxo(tab = tab, # aggregate my data at the targeted tax level
                                     tax.lvl = tax.lvl,
                                     tax.table = SDPs_16s@tax_table) %>%
  apply(1,sum) %>% # sum the number of reads
  sort %>% # sort in ascending order
  identity
top_class <- class_tot%>%
  tail(n = 9)%>%
  reshape2::melt()%>%
  mutate(class=rownames(.))%>%
  mutate(class=ifelse(nchar(class)<1,"Unidentified",class))

# run pcoa=======================================

dtmp <- file2meco::phyloseq2meco(SDPs_16s)
dtmp$tax_table%<>%
  mutate(prev=rowSums(dtmp$otu_table>0))

dtmp$otu_table <- dtmp$otu_table[rowSums(dtmp$otu_table)>1,]#remove singletons
dtmp$otu_table <- t(labdsv::hellinger(t(dtmp$otu_table)))  #perform hellinger transformation
orditmp <- phyloseq::ordinate(file2meco::meco2phyloseq(dtmp), "PCoA", "horn") # run pcoa with Morisita-Horn distance

# get relative eigenvalues (% variance associated to different axes)
expl1 <- orditmp$values$Relative_eig[1] #axis1
expl2 <- orditmp$values$Relative_eig[2] #axis2
expl3 <- orditmp$values$Relative_eig[3] #axis3

###Samples============================================
# Set climate zones color
climate_palette <- paletteer::paletteer_d("MetBrewer::Austria")[1:3]
names(climate_palette) <- c("Continental hot summer","Continental cold/Tundra","Temperate")


#check the hypothesis that dispersion is roughly the same between my groups in multivariate space
dis <- vegan::vegdist(labdsv::hellinger(t(dtmp$otu_table)),method="horn") #we create a distance matrix between our individuals
mod <- vegan::betadisper(dis,dtmp$sample_table$koppen_clim) #we calculate the multivariate dispersions i.e. average distance to centroids
p_dist2grpcentroid <- 
  mod$distances%>%
  reshape2::melt()%>%
  mutate(Sample=rownames(.))%>%
  left_join(dtmp$sample_table)%>%
  ggplot(aes(fill=koppen_clim,x=value))+
  geom_density(alpha=.5)+
  xlab("Distance to group centroïd")+
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = climate_palette,name="Climatic zone")+
  theme(text = element_text(face = "bold"))
p_dist2grpcentroid


# We see that distance to group centroïd is roughly the same so that it's ok to perform a permanova to test for shifts in group centroids between climatic zones.

##run permanova===========
adon.res1 <- vegan::adonis2(labdsv::hellinger(t(dtmp$otu_table))~dtmp$sample_table$koppen_clim,
                            perm = 1000, na.rm = T, method = "horn",p.adjust.methods="holm" ) 


adon.res1%>% #save # check path or remove
  flextable::flextable()%>%
  flextable::save_as_docx(path=paste0("figures_and_results/tables/permanova_comunities.docx"))

pair.adon.res1 <- pairwiseAdonis::pairwise.adonis(labdsv::hellinger(t(dtmp$otu_table)), #run pairwise permanova
                                                  dtmp$sample_table$koppen_clim,
                                                  p.adjust.m = "holm",
                                                  sim.method = "horn", perm=1000)
pair.adon.res1%>% #save
  dplyr::select(-sig)%>%
  flextable::flextable(cwidth = c(1,.5,1,1,1,1,1))%>%
  flextable::save_as_docx(path=paste0("figures_and_results/tables/permanova_comunities_pairwise.docx"))


# Extract sample coordinates on 3 first dimensions
dim1 <- orditmp$vectors[,1]
dim2 <- orditmp$vectors[,2]
dim3 <- orditmp$vectors[,3]
Sample <- names(orditmp$vectors[,1])
dtmp$sample_table%<>%
  mutate(dim1=dim1,
         dim2=dim2,
         dim3=dim3)
# create df and merge with sample info

df_ordihorn <- dtmp$sample_table
df_ordihorn%<>%
  group_by(koppen_clim)%>%
  mutate(medx=median(dim1), #median of the coordinates
         medy=median(dim2),
         medy3=median(dim3))
# plot Sample PCoA
pcoa_horn <- ggplot(df_ordihorn,aes(x=dim1,y=dim2,color=koppen_clim))+
  stat_centseg()+
  # add point corresponding to class median position
  geom_point(mapping=aes(medx,medy,
                         color=koppen_clim), #bigger size for top 10 classes
             alpha=.5,shape=21,fill="white",key_glyph=draw_key_point, stroke = 1,size=3)+
  ylab(paste0("Axis 2 (",round(expl2,digits = 3)*100,"%)"))+
  xlab(paste0("Axis 1 (",round(expl1,digits = 3)*100,"%)"))+
  # ggConvexHull::geom_convexhull(aes(fill = koppen_clim, color = koppen_clim),alpha=0.2)+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  scale_color_manual(values = climate_palette,name="Climatic zone")+
  theme_bw()+
  theme(text=element_text(face="bold"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(fill=climate_palette[c(2,1,3)],alpha=1)))


pcoa_horn13 <- ggplot(df_ordihorn,aes(x=dim1,y=dim3,color=koppen_clim))+
  stat_centseg()+
  # add point corresponding to class median position
  geom_point(mapping=aes(medx,medy3,
                         color=koppen_clim), #bigger size for top 10 classes
             alpha=.5,shape=21,fill="white", stroke = 1,size=3)+
  theme_bw()+
  theme(text=element_text(face="bold"),
        panel.grid = element_blank())+
  ylab(paste0("Axis 3 (",round(expl3,digits = 3)*100,"%)"))+
  xlab(paste0("Axis 1 (",round(expl1,digits = 3)*100,"%)"))+
  # ggConvexHull::geom_convexhull(aes(fill = koppen_clim, color = koppen_clim),alpha=0.2)+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  scale_color_manual(values = climate_palette,name="Climatic zone")+
  theme(legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(fill=climate_palette[c(2,1,3)],alpha=1)))




### Taxa========================

#get taxa coordinates for axes 1,2,3
pordi_asv <-  phyloseq::plot_ordination(file2meco::meco2phyloseq(dtmp), orditmp,type = "taxa",axes = c(1,2))
pordi_asv2 <-  phyloseq::plot_ordination(file2meco::meco2phyloseq(dtmp), orditmp,type = "taxa",axes = c(1,3))

#create dataframe
df <- pordi_asv$data%>%
  mutate(nb_reads=rowSums(dtmp$otu_table))%>% #compute number of reads
  arrange(prev)%>% #order by prevalence
  mutate(ASV_id=forcats::fct_inorder(ASV_id))%>% #level factor
  left_join(select(pordi_asv2$data,c(2,3)))%>% #add axis 3
  group_by(Class)%>% # groupby class
  mutate(maxprev=max(prev), # prevalence of the most prevalent ASV
         propreads=sum(nb_reads)/sum(dtmp$otu_table[rowSums(dtmp$otu_table)>1,]), # proportion of total reads
         medx=median(Axis.1), #median of the coordinates
         medy=median(Axis.2),
         medy3=median(Axis.3))%>%
  ungroup()%>%
  mutate(Class=gsub("c__","",Class)) %>%
  mutate(class2plot=ifelse(Class%in%top_class$class[-c(which(top_class$class=="Unidentified"))],Class,"Others"))%>% # Keep only class that we want and tag the rest as "Others"
  mutate(is2plot=ifelse(Class%in%top_class$class[-c(which(top_class$class=="Unidentified"))],'yup',"nope"))%>% # Create variable saying if to plot or "others"
  mutate(class2plot=forcats::fct_relevel(class2plot,"Others"))%>% #place Others as first factor
  arrange(class2plot) # arrange so that the dots are drawn in right order


dperm_taxa <- df%>%
  select(Axis.1,Axis.2,Axis.3,class2plot)%>%
  filter(class2plot!="Others")
  
#check the hypothesis that dispersion is roughly the same between my groups in multivariate space
dis <- vegan::vegdist(dperm_taxa[,c(1:3)],method="euclidean") #we create a distance matrix between our individuals
mod <- vegan::betadisper(dis,dperm_taxa$class2plot) #we calculate the multivariate dispersions i.e. average distance to centroids

p_dist2grpcentroid <- 
  mod$distances%>%
  reshape2::melt()%>%
  mutate(ASV_id=filter(df,class2plot!="Others")$ASV_id)%>%
  left_join(df)%>%
  ggplot(aes(fill=class2plot,x=value,color=class2plot))+
  geom_density(alpha=.4)+
  xlab("Distance to group centroïd")+
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c(paletteer::paletteer_d("MoMAColors::Klein")[1:8]),name="Class")+
  scale_color_manual(values=c(paletteer::paletteer_d("MoMAColors::Klein")[1:8]),name="Class")+
  theme(text = element_text(face = "bold"))+
  guides(fill=guide_legend(override.aes = list(alpha=1)))

p_dist2grpcentroid

# We see that distance to group centroïd is roughly homogeneous so that it's ok to perform a permanova to test for shifts in group centroids between climatic zones.
#run permanova
adon.res1 <- vegan::adonis2(dperm_taxa[,c(1:3)]~dperm_taxa$class2plot,
                            perm = 1000, na.rm = T, method = "euclidean",p.adjust.methods="holm" ) 
adon.res1%>% #save
  flextable::flextable()%>%
  flextable::save_as_docx(path=paste0("figures_and_results/tables/permanova_taxa.docx"))

pair.adon.res1 <- pairwiseAdonis::pairwise.adonis(dperm_taxa[,c(1:3)], #run pairwise permanova
                                                  dperm_taxa$class2plot,
                                                  p.adjust.m = "holm",
                                                  sim.method = "euclidean", perm=1000)
pair.adon.res1%>% #save
  dplyr::select(-sig)%>%
  flextable::flextable(cwidth = c(1,.5,1,1,1,1,1))%>%
  flextable::save_as_docx(path=paste0("figures_and_results/tables/permanova_taxa_pairwise.docx"))
##Plot=====================

# Axis bxplot
axis1 <- df%>%
  filter(class2plot!="Others")%>% # remove "others" classes as grouping them makes no sens
  
  ggplot(aes(y=class2plot,x=Axis.1,fill=class2plot))+ #plot coordinates according to classes
  geom_violin(key_glyph=draw_key_dotplot,color=NA)+
  geom_boxplot(width=0.2,fill="white",color="black",outliers = F)+
  scale_fill_manual(values=c(paletteer::paletteer_d("MoMAColors::Klein")[1:8]),name="Class",guide="none")+
  #set lims exactly as those of pcoa
  xlim(layer_scales(pcoa_horn)$x$range$range)+
  #set theme
  theme_bw()+
  theme(text=element_text(face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  geom_vline(xintercept=0,lty=2)


axis2 <- df%>%
  filter(class2plot!="Others")%>%
  ggplot(aes(class2plot,Axis.2,fill=class2plot))+
  geom_violin(key_glyph=draw_key_dotplot,color=NA)+
  geom_boxplot(width=0.2,fill="white",color="black",outliers = F)+
  theme_bw()+
  scale_fill_manual(values=paletteer::paletteer_d("MoMAColors::Klein")[1:8],name="Class",guide="none")+
  ylim(layer_scales(pcoa_horn)$y$range$range)+
  theme(text=element_text(face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  geom_hline(yintercept=0,lty=2)

axis3 <- df%>%
  filter(class2plot!="Others")%>%
  ggplot(aes(class2plot,Axis.3,fill=class2plot))+
  geom_violin(key_glyph=draw_key_dotplot,color=NA)+
  geom_boxplot(width=0.2,fill="white",color="black",outliers = F)+
  theme_bw()+
  ylim(layer_scales(pcoa_horn13)$y$range$range)+
  scale_fill_manual(values=paletteer::paletteer_d("MoMAColors::Klein")[1:8],name="Class",guide="none")+
  theme(text=element_text(face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  geom_hline(yintercept=0,lty=2)

#compose plot
pcoataxfull <- (axis1+ggpubr::get_legend(axis1+
                                           guides(fill=guide_legend(reverse = T)))+pcoa_horn+axis2+pcoa_horn13+axis3)+
  #set layout
  plot_layout(guides="collect",heights=c(3,4,4),widths = c(4,3))&
  # set theme
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "bottom")

pcoataxfull
