#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 3A - Enzymatic biogeographic patterns and drivers
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Import data
enz_dat <- read.csv('enzyme_data.csv', row.names=1)


###
#### Response to Koppen climate ----
###

#BG
mod <- aov(log(enz_dat$BG+1) ~enz_dat$koppen_clim)
hist(mod$residuals)
summary(mod, type ='III') #p = 0.0123*
boxplot(enz_dat$BG ~enz_dat$koppen_clim)
TukeyHSD(mod)
stat_bg <- HSD.test(mod, "enz_dat$koppen_clim", group=TRUE)
stat_bg$groups$groups

#NAG
mod <- aov(log(enz_dat$NAG+1) ~enz_dat$koppen_clim)
hist(mod$residuals)
summary(mod, type ='III') #p < 0.001***
boxplot(log(enz_dat$NAG+1) ~enz_dat$koppen_clim)
TukeyHSD(mod)
stat_nag <- HSD.test(mod, "enz_dat$koppen_clim", group=TRUE)
stat_nag

#LAP
mod <- aov(log(enz_dat$LAP+1) ~enz_dat$koppen_clim)
hist(mod$residuals)
summary(mod, type ='III') #p = 0.002*
boxplot(log(enz_dat$LAP+1) ~enz_dat$koppen_clim)
TukeyHSD(mod)
stat_lap <- HSD.test(mod, "enz_dat$koppen_clim", group=TRUE)
stat_lap$groups$groups

#AP
mod <- aov(log(enz_dat$AP+1) ~enz_dat$koppen_clim)
hist(mod$residuals)
summary(mod, type ='III') #p < 0.001***
boxplot(log(enz_dat$AP+1) ~enz_dat$koppen_clim)
TukeyHSD(mod)
stat_ap <- HSD.test(mod, "enz_dat$koppen_clim", group=TRUE)
stat_ap


## Boxplots

enz_bg <- ggplot(enz_dat, aes(x=koppen_clim,y=BG))+
  geom_boxplot(outlier.shape = NA,aes(fill=koppen_clim))+
  geom_jitter(aes(col=koppen_clim), width = .1, height = 0, alpha = 0.4)+
  geom_text(data= as.data.frame(stat_bg$groups), aes(x = row.names(stat_bg$groups), y = 1500, label = stat_bg$groups$groups))+
  scale_fill_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  scale_color_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  theme_bw()+
  ylab(bquote(bold("BG activity (nmol.h"^-1~".g"[DM]~")")))+
  theme(text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
enz_nag <- ggplot(enz_dat, aes(x=koppen_clim,y=NAG))+
  geom_boxplot(outlier.shape = NA,aes(fill=koppen_clim))+
  geom_jitter(aes(col=koppen_clim), width = .1, height = 0, alpha = 0.4)+
  geom_text(data= as.data.frame(stat_nag$groups), aes(x = row.names(stat_nag$groups), y = 2000, label = stat_nag$groups$groups))+
  scale_fill_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  scale_color_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  theme_bw()+
  ylab(bquote(bold("NAG activity (nmol.h"^-1~".g"[DM]~")")))+
  theme(text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
enz_lap <- ggplot(enz_dat, aes(x=koppen_clim,y=LAP))+
  geom_boxplot(outlier.shape = NA,aes(fill=koppen_clim))+
  geom_jitter(aes(col=koppen_clim), width = .1, height = 0, alpha = 0.4)+
  geom_text(data= as.data.frame(stat_lap$groups), aes(x = row.names(stat_lap$groups), y = 200, label = stat_lap$groups$groups))+
  scale_fill_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  scale_color_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  theme_bw()+
  ylab(bquote(bold("LAP activity (nmol.h"^-1~".g"[DM]~")")))+
  theme(text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
enz_ap <- ggplot(enz_dat, aes(x=koppen_clim,y=AP))+
  geom_boxplot(outlier.shape = NA,aes(fill=koppen_clim))+
  geom_jitter(aes(col=koppen_clim), width = .1, height = 0, alpha = 0.4)+
  geom_text(data= as.data.frame(stat_ap$groups), aes(x = row.names(stat_ap$groups), y = 7000, label = stat_ap$groups$groups))+
  scale_fill_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  scale_color_manual(values=c("#7C8AB5",  "#C45858","#6EB689"))+
  theme_bw()+
  ylab(bquote(bold("AP activity (nmol.h"^-1~".g"[DM]~")")))+
  theme(text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
gridExtra::grid.arrange(enz_bg, enz_nag, enz_lap, enz_ap)

###
#### Environmental drivers ----
###

## Import data
enz_env <- read.csv("Enzymes_and_Env_drivers.csv", row.names=1) #input file

### Importance plot

# Load data
df <- enz_env[,-c(1,3,4)] #reproduce this step foir each enzyme by removing the other enzymes. 

# Calculate Spearman rho for each predictor vs response
response_var <- "NAG"  # replace with your actual response variable
predictors <- setdiff(names(df), response_var)

# Initialize result storage
results <- data.frame(
  Variable = predictors,
  SpearmanRho = NA,
  P_value = NA
)

for (v in predictors) {
  test <- cor.test(df[[v]], df[[response_var]], method = "spearman")
  results$SpearmanRho[results$Variable == v] <- test$estimate
  results$P_value[results$Variable == v] <- test$p.value
}

# Add derived columns
results$AbsRho <- abs(results$SpearmanRho)
results$Significance <- cut(results$P_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", ""))
results$SigGroup <- ifelse(results$P_value < 0.05, "Significant", "Non-significant")

# Sort by importance
results <- results[order(-results$AbsRho), ]

# Plot

imp_bg <- ggplot(results, aes(x = reorder(Variable, AbsRho), y = AbsRho)) +
  geom_bar(data = subset(results, SigGroup == "Non-significant"),
           stat = "identity", fill = "lightgray") +
  geom_bar(data = subset(results, SigGroup == "Significant"),
           aes(fill = AbsRho), stat = "identity") +
  geom_text(aes(label = Significance), hjust = -0.2, size = 5) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  coord_flip() +
  labs(title = "BG",
       y = expression(paste("|Spearman ", rho, "|")),
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

imp_nag <- ggplot(results, aes(x = reorder(Variable, AbsRho), y = AbsRho)) +
  geom_bar(data = subset(results, SigGroup == "Non-significant"),
           stat = "identity", fill = "lightgray") +
  geom_bar(data = subset(results, SigGroup == "Significant"),
           aes(fill = AbsRho), stat = "identity") +
  geom_text(aes(label = Significance), hjust = -0.2, size = 5) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  coord_flip() +
  labs(title = "NAG",
       y = expression(paste("|Spearman ", rho, "|")),
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

imp_lap <- ggplot(results, aes(x = reorder(Variable, AbsRho), y = AbsRho)) +
  geom_bar(data = subset(results, SigGroup == "Non-significant"),
           stat = "identity", fill = "lightgray") +
  geom_bar(data = subset(results, SigGroup == "Significant"),
           aes(fill = AbsRho), stat = "identity") +
  geom_text(aes(label = Significance), hjust = -0.2, size = 5) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  coord_flip() +
  labs(title = "LAP",
       y = expression(paste("|Spearman ", rho, "|")),
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

imp_ap <- ggplot(results, aes(x = reorder(Variable, AbsRho), y = AbsRho)) +
  geom_bar(data = subset(results, SigGroup == "Non-significant"),
           stat = "identity", fill = "lightgray") +
  geom_bar(data = subset(results, SigGroup == "Significant"),
           aes(fill = AbsRho), stat = "identity") +
  geom_text(aes(label = Significance), hjust = -0.2, size = 5) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  coord_flip() +
  labs(title = "AP",
       y = expression(paste("|Spearman ", rho, "|")),
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

gridExtra::grid.arrange(imp_bg, imp_nag, imp_lap, imp_ap)

