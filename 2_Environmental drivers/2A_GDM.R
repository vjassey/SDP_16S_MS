#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 2A - General dissimilarity modelling
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(gdm)
library(ggplot2)
library(viridis)
library(patchwork)

## Import data
load("SDPs_16s.RData") # input file
asvs <- as.data.frame(as.matrix(t(SDPs_16s@otu_table)))
env_asv <- read.csv('env_predictors_cleaned.csv', row.names=1)
coordinates <- read.csv('coordinates_cleaned.csv', row.names=1)

# Transform data
env_asv_logged <- env_asv

## Apply log transformation
env_asv_logged[cols_to_log_transform] <- lapply(env_asv_logged[cols_to_log_transform], function(x) log(x+1))
env_asv_logged$koppen_clim <- coordinates$koppen_clim


#select env variables
DataSet <- env_asv_logged %>% dplyr::select(swe_lt , tmmn_lt, height, pr_lt , ph, gpp_lt, nonTreeCover, soil_moisture_lt)
DataSet <- decostand(DataSet, 'standardize')
DataSet <- cbind(coordinates[,c(2,3)], DataSet)

# Site by species matrix ('pa' data)
sp.pa <- decostand(bio.b, "pa")
range(colSums(sp.pa))
sp_sdm <- cbind(sp.pa, env[,c(2,3,4)])

# Site by species matrix ('ab' data) + coordinates
sp_sdm <- cbind(asvs, coordinates[,c(2,3)])

#Remove empty rows (community with only = 0)
#sp_sdm <- sp_sdm[rowSums(sp_sdm[,1:1049])>0,]

#Site Column
sp_sdm$Sample <- env_asv$Sample
range(colSums(asvs))

# get columns with xy, site ID, and species data (abundance data)
sppTab <- sp_sdm[,1:15218]
sppTab <- sppTab %>% select_if(colSums(.) > 0) # remove empty columns
sppTab <- sppTab[rowSums(sppTab[])>0,] #remove empty rows
Sample <- env_asv$Sample
sppTab <- cbind(coordinates[,c(2, 3)], Sample, sppTab[,1:15218]) # get columns with env. data and xy-coordinates Site, Lon and Lat; in this order!

#pa
sppTab <- cbind(env[,c(3, 4)], Sample, sppTab[,1:5914]) # get columns with env. data and xy-coordinates Site, Lon and Lat; in this order!


#envTab
envTab <- cbind(Sample, DataSet)
envTab$X <- as.numeric(envTab$X); envTab$Y <- as.numeric(envTab$Y)
str(envTab)

# Compute geographic distance matrix using geodesic method to account for earth's curvature (GDM is not)
mat_geo_gdm <- geodist::geodist(coordinates[,c("X","Y")],measure = 'geodesic')/1000 # /1000 for km distance
rownames(mat_geo_gdm) <- colnames(mat_geo_gdm) <- env_asv$Sample
mat_geo_gdm <- cbind("Sample"=rownames(mat_geo_gdm),mat_geo_gdm)

# x-y species list example
gdmTab<- formatsitepair(sppTab, dist="bray", abundance=F, bioFormat=1, siteColumn = "Sample", XColumn="X", YColumn="Y", 
                        sppColumn=NA, predData=envTab, abundColumn=NA, distPreds =NULL) 
gdm.1 <- gdm(gdmTab, geo=T) 
plot(gdm.1)
par(ask=F)

# Variable importance using permutations
varimp = gdm.varImp(gdmTab, geo=T, nPerm=10, parallel = T, cores=5) # try several perm and uses violin in the plot?
varimpImp <- varimp$`Predictor Importance`
varimpImp$preds <- rownames(varimpImp)
varimpImp <- arrange(varimpImp, desc(`All predictors`)) #arrange by decreasing importance
colnames(varimpImp) <- c('importance', 'variable') 
#par(ask = FALSE) # to remove Hit <Return> to see next plot:

###
### plot variance partitioning
###

# Percent deviance explained by the full model =  26.979 

#The relative importance of climatic predictors was estimated as the sum of I-spline basis functions 
# Maximum height of each response curve (Fitzpatrick and Keller, 2015; Robroek et al. 2017).
summary(gdm.1)
gdm.1.splineDat <- isplineExtract(gdm.1)

#save splines
max_splines <-as.data.frame(as.matrix(gdm.1.splineDat$y))

# Maximum spline values for each column
max_val <- as.data.frame(apply(max_splines,2,max))
colnames(max_val) <- 'max' 
max_val$variable <- row.names(max_val)
max_val[10,] <- c(sum(max_val$max[c(2,3,5)]), "climate")
max_val[11,] <- c(sum(as.numeric(max_val$max[c(4,6,7,8)])), "local")

# Partition
part_gdm <- max_val[c(1,10,11),]
part_gdm$max <- as.numeric(part_gdm$max)
part_gdm$type <- 'gdm_ab'

#barplot of variance partitioning
varpart_p <- ggplot(part_gdm, aes(fill=variable, y=max, x=type)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + scale_fill_viridis(discrete = T)

#' 
#' 
#' Plot GDM splines
## ----------------------------------------------------------------------------------------------------------------------------------------------------
# get predictor significance
gdm_pred_signif <- varimp$`Predictor p-values`%>%
  mutate(Var2=gsub("_[a-z]{2}","",rownames(.)))%>%
  rename(pval=`All predictors`)%>%
  mutate(issignif=ifelse(pval<0.05,"yup","nope"))

# get spline data
gdm_spline_data <- gdm::isplineExtract(gdm.1)

# format data for plot
dgdm_var <- gdm_spline_data$x%>% # get variable variation range
  reshape2::melt()%>%
  rename(valuex=value)%>%
  #get variable with non null GDM effect
  left_join(reshape2::melt(gdm_spline_data$y[,which(colSums(gdm_spline_data$y)!=0)]))%>% 
  mutate(Var2=gsub("_[a-z]{2}","",Var2))%>%
  rename(valuey=value)%>%
  #left_join(names_df_clean)%>%  #join clean names
  #mutate(midname=ifelse(is.na(midname),"Geographic distance",midname))%>% #rename
  #mutate(midname=forcats::fct_relevel(midname,unique(names_df_clean$midname)))%>% #relevel
  filter(Var2%in%gsub("_[a-z]{2}","",names(which(colSums(gdm_spline_data$y)!=0))))%>%
  left_join(gdm_pred_signif)

#plot variable spline
gdm_plot_var <- dgdm_var %>%
  group_by(Var2)%>% # group by variable
  #normalize variable unit between 0-1
  mutate(valuexn=((valuex - min(valuex, na.rm = T))/(max(valuex, na.rm = T) - min(valuex, na.rm = T))))%>% 
  mutate(mxy=max(valuey))%>% # get maximal value
  ungroup()%>%
  arrange(desc(mxy))%>% # arrange by maximal value
  mutate(Var2=forcats::fct_inorder(Var2))%>% # order levels by max value
  ggplot(aes(x=valuexn, y=valuey,color=Var2,lty=issignif))+ #plot lines
  geom_line(linewidth=1.5)+
  theme_bw()+
  scale_color_manual(values=paletteer::paletteer_d("rcartocolor::TealRose"))+
  scale_linetype_manual(values=c(2,1),guide="none")+
  guides(color=guide_legend(ncol=3,position = "inside"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  scale_x_continuous(expand = c(0,0))+
  xlab("0-1 scaled variables")+
  ylab("GDM effect")+
  theme(legend.position.inside = c(.5,.9),
        axis.text.x = element_text(angle=45,hjust = 1),
        text=element_text(face="bold"),
        legend.background = element_rect(color='black'),
        legend.title = element_blank())


varpart_p + gdm_plot_var 

#' 
#' GDM Geographic only ----
## ----------------------------------------------------------------------------------------------------------------------------------------------------

# GDM AB versus Geographic distance only
gdmTab_geo<- formatsitepair(sppTab, dist="bray", abundance=F, bioFormat=1, siteColumn = "Sample", XColumn="X", YColumn="Y", 
                            sppColumn=NA, predData=envTab[,c(1:3)], abundColumn=NA, distPreds =NULL) 
#for plotting within/between
envTab2 <- cbind(envTab, coordinates$koppen_clim)
gdmTab_geo2<- formatsitepair(sppTab, dist="bray", abundance=F, bioFormat=1, siteColumn = "Sample", XColumn="X", YColumn="Y", 
                             sppColumn=NA, predData=envTab2[,c(1:3, 12)], abundColumn=NA, distPreds =NULL) 
gdm.ab_geo <- gdm(gdmTab_geo, geo=T) 
plot(gdm.ab_geo)

# GDM P/A versus Geographic distance only
sppTab_pa <- decostand(sp_sdm[,1:15218], 'pa')
sppTab_pa <- sppTab_pa %>% select_if(colSums(.) > 1) # remove empty columns
sppTab_pa <- sppTab_pa[rowSums(sppTab_pa[])>0,] #remove empty rows
sppTab_pa <- cbind(coordinates[,c(2, 3)], Sample, sppTab_pa)

gdmTab_geo_pa<- formatsitepair(sppTab_pa, dist="jaccard", abundance=F, bioFormat=1, siteColumn = "Sample", XColumn="X", YColumn="Y", 
                               sppColumn=NA, predData=envTab[,c(1:3)], abundColumn=NA, distPreds =NULL) 
gdm.pa_geo <- gdm(gdmTab_geo_pa, geo=T) 
plot(gdm.pa_geo)


# 
# 
# Plot observation diss versus Geographic distance of the GDM ----
## ----------------------------------------------------------------------------------------------------------------------------------------------------
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

### AB GDM versus geographic distance
#get data
dgdm_pred <- data.frame(pred=gdm.ab_geo$predicted,
                        obs=gdm.ab_geo$observed,
                        eco=gdm.ab_geo$ecological)


#get density for plotting
dgdm_pred$density_ecoobs <- get_density(dgdm_pred$obs, dgdm_pred$eco, n = 100)
dgdm_pred$density_obspred <- get_density(dgdm_pred$obs, dgdm_pred$pred, n = 100)

#get exp curve to plot over data
overlayX <- seq(from = min(dgdm_pred$eco), 
                to = max(dgdm_pred$eco), 
                length = 171)
overlayY <- 1 - exp(-overlayX)
df_over <- data.frame(x=overlayX,y=overlayY)

#get koppen_clim info
dgdm_pred$s1.koppen_clim <- gdmTab_geo2$`s1.coordinates$koppen_clim`
dgdm_pred$s2.koppen_clim <- gdmTab_geo2$`s2.coordinates$koppen_clim`

# Add area info
dgdm_pred <- dgdm_pred %>%
  mutate(
    Area1 = as.character(s1.koppen_clim),
    Area2 = as.character(s2.koppen_clim),
    ComparisonType = ifelse(Area1 == Area2, "Within", "Between")
  )

# Plot prediction vs obs with Within/Between diss
gdm_plot_obs <- ggplot(dgdm_pred, aes(x = eco, y = obs, color = ComparisonType)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Geographic Distance",
    y = "Observed Dissimilarity",
    color = "Bioclimatic Comparison"
  ) +
  scale_color_brewer(palette = "Set1")+
  geom_line(data=df_over,mapping=aes(x=x,y=y),color="#BCB584",lwd=1.5)+
  ylim(c(0.5,1))+ 
  xlab("Geographic distance")+
  ylab("Observed dissimilarity")+
  ggpubr::theme_classic2()+
  theme(legend.position = c(.9,.25),
        legend.background = element_rect(color="black"),
        legend.title = element_text( size=7),
        legend.text=element_text(size=7),
        legend.key.height = unit(8, 'pt'),
        legend.key.width = unit(4,"pt"))

# Plot prediction vs obs
gdm_plot_obs <-
  #ecological distance vs observed biodistance
  ggplot()+
  geom_point(dgdm_pred,mapping=aes(x=eco,
                                   y=obs,
                                   color=density_ecoobs),size=1)+
  geom_line(data=df_over,mapping=aes(x=x,y=y),color="#BCB584",lwd=1.5)+
  ylim(c(0.5,1))+ 
  scale_color_gradientn(colors = paletteer::paletteer_d("fishualize::Acanthostracion_polygonius"),
                        name="Density")+
  xlab("Geographic distance")+
  ylab("Observed dissimilarity")+
  ggpubr::theme_classic2()+
  theme(legend.position = c(.9,.25),
        legend.background = element_rect(color="black"),
        legend.title = element_text( size=7),
        legend.text=element_text(size=7),
        legend.key.height = unit(8, 'pt'),
        legend.key.width = unit(4,"pt"))
gdm_plot_obs


### P/A GDM versus geographic distance
#get data
dgdm_pred_pa <- data.frame(pred=gdm.pa_geo$predicted,
                           obs=gdm.pa_geo$observed,
                           eco=gdm.pa_geo$ecological)

#get density for plotting
dgdm_pred_pa$density_ecoobs <- get_density(dgdm_pred_pa$obs, dgdm_pred_pa$eco, n = 100)
dgdm_pred_pa$density_obspred <- get_density(dgdm_pred_pa$obs, dgdm_pred_pa$pred, n = 100)

#get exp curve to plot over data
overlayX_pa <- seq(from = min(dgdm_pred_pa$eco), 
                   to = max(dgdm_pred_pa$eco), 
                   length = 171)
overlayY_pa <- 1 - exp(-overlayX_pa)
df_over_pa <- data.frame(x=overlayX_pa,y=overlayY_pa)

#get koppen_clim info
dgdm_pred_pa$s1.koppen_clim <- gdmTab_geo2$`s1.coordinates$koppen_clim`
dgdm_pred_pa$s2.koppen_clim <- gdmTab_geo2$`s2.coordinates$koppen_clim`

# Add area info
dgdm_pred_pa <- dgdm_pred_pa %>%
  mutate(
    Area1 = as.character(s1.koppen_clim),
    Area2 = as.character(s2.koppen_clim),
    ComparisonType = ifelse(Area1 == Area2, "Within", "Between")
  )

# Plot prediction vs obs with Within/Between diss
gdm_plot_obs_pa <- ggplot(dgdm_pred_pa, aes(x = eco, y = obs, color = ComparisonType)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Geographic Distance",
    y = "Observed Dissimilarity",
    color = "Bioclimatic Comparison"
  ) +
  scale_color_brewer(palette = "Set1")+
  geom_line(data=df_over_pa,mapping=aes(x=x,y=y),color="#BCB584",lwd=1.5)+
  ylim(c(0.5,1))+
  xlab("Geographic distance")+
  ylab("Observed dissimilarity")+
  ggpubr::theme_classic2()+
  theme(legend.position = c(.9,.25),
        legend.background = element_rect(color="black"),
        legend.title = element_text( size=7),
        legend.text=element_text(size=7),
        legend.key.height = unit(8, 'pt'),
        legend.key.width = unit(4,"pt"))
gdm_plot_obs_pa

# Plot prediction vs obs
gdm_plot_obs_pa <-
  #ecological distance vs observed biodistance
  ggplot()+
  geom_point(dgdm_pred_pa,mapping=aes(x=eco,
                                      y=obs,
                                      color=density_ecoobs),size=1)+
  geom_line(data=df_over_pa,mapping=aes(x=x,y=y),color="#BCB584",lwd=1.5)+
  ylim(c(0.5,1))+ 
  scale_color_gradientn(colors = paletteer::paletteer_d("fishualize::Acanthostracion_polygonius"),
                        name="Density")+
  xlab("Geographic distance")+
  ylab("Observed dissimilarity")+
  ggpubr::theme_classic2()+
  theme(legend.position = c(.9,.25),
        legend.background = element_rect(color="black"),
        legend.title = element_text( size=7),
        legend.text=element_text(size=7),
        legend.key.height = unit(8, 'pt'),
        legend.key.width = unit(4,"pt"))
gdm_plot_obs_pa  

#observed biodistance vs pred
ggplot(dgdm_pred,
       aes(x=pred,
           y=obs,
           color=density_obspred))+
  geom_point(size=1)+
  geom_abline (slope=1, color="#BCB584",lwd=1.5)+
  scale_color_gradientn(colors = paletteer::paletteer_d("fishualize::Acanthostracion_polygonius"),
                        name="Density")+
  xlab("Predicted dissimilarity")+
  ylab("Observed dissimilarity")+
  annotate(geom = "label",x=.5,y=1,
           label=paste0("Deviance explained: ",round(gdm_fitted$explained,digits=2),"%\nAvg obs dissimilarity: ",round(gdm_fitted$intercept,digits=2)),
           hjust=0,fontface="bold",vjust=1,
           fill=fill_alpha("lightgrey", .8),
           color="black")+
  # ylim(c(.2,1))+
  # xlim(c(.2,1))+
  ggpubr::theme_classic2()+
  theme(legend.position = c(.9,.25),
        legend.background = element_rect(color="black"),
        legend.title = element_text( size=7), 
        legend.text=element_text(size=7),
        legend.key.height = unit(8, 'pt'),
        legend.key.width = unit(4,"pt"))&
  theme(text=element_text(face="bold"))

gridExtra::grid.arrange(gdm_plot_obs, gdm_plot_obs_pa)  




