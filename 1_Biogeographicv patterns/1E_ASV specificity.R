#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 1E - ASV specificity
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## This code is adapted from Ezzat, L., Peter, H., Bourquin, M. et al. 
#Diversity and biogeography of the bacterial microbiome in glacier-fed streams. 
#Nature 637, 622–630 (2025). https://doi.org/10.1038/s41586-024-08313-z

library(speedyseq)
library(phyloseq)
library(tidyverse)
library(phyloseqCompanion)
library(strider)

## Load ASV and taxonomy tables
load("SDPs_16s.RData") # input file

## Create phyloseq object
merged_NOMIS_DEBLUR <- SDPs_16s
metadata_nomis_inext <- sample.data.frame(merged_NOMIS_DEBLUR)
prune_asv16s_df <- as.matrix(SDPs_16s@otu_table)
row_sums <- rowSums(prune_asv16s_df )

## Subset the dataframe to remove rows where the sum is zero
prune_asv16s_df <- prune_asv16s_df[row_sums > 0, ] ### 15213 ASVs and 171 SDPs
metadata_nomis <- sample.data.frame(merged_NOMIS_DEBLUR)
asv_df <- as.data.frame(t(otu_table(prune_asv16s_df, taxa_are_rows=T)))

## Here we are investigating the prevalence of ASV across the dataset
howmanyasv<- as.data.frame(colSums(asv_df != 0))
colnames(howmanyasv) <- c("Count_nb")
howmanyasv$ASV <- rownames(howmanyasv)
rownames(howmanyasv) <- NULL

## Melt asv table
asvdfmelt <- reshape2::melt(as.matrix(asv_df))
## Keep only values that are > 0
asvdfmelt <- asvdfmelt[asvdfmelt$value >0,]

##### Specific ASVs
#up + down - List of ASVs occurring in only one SPD: 9304, 9304/15218=61.1%
endemic_oneSPD <- howmanyasv %>% 
  filter(Count_nb == 1) 
#write.csv(endemic_oneSPD, "MAPP_asv16S_unique_list.csv")

endemic_oneSPD_n <- howmanyasv %>% 
  filter(Count_nb == 1) %>%
  summarise(n = n())

## total number of specific ASVs found in one single SPD
endemic_oneSPD_n/length(unique(row.names(prune_asv16s_df)))#61.1%

# value is the nb of count and count_nb nb of samples 
merge_equalone <- merge(endemic_oneSPD, asvdfmelt, by.x="ASV",by.y="Var2")
## assign the name of the mountain range
merge_equalone$MR <- vapply((merge_equalone$Var1), function(x) metadata_nomis$koppen_clim[metadata_nomis$Sample == x], FUN.VALUE = character(1))
## rename column
colnames(merge_equalone) <- c("ASV","prev","glname","nb_count","MR")

## Nb of ASVs that are range-specific 
endemic_MR_prev <- merge_equalone %>% group_by(MR) %>% summarize(prev=sum(prev))
endemic_oneSPD_plot <- merge_equalone %>% group_by(MR) 
endemic_oneSPD_plot<- endemic_oneSPD_plot[c("MR","ASV","prev")]
endemic_oneSPD_plot$Color <- "C_uniquetoone"

## Calculate unique ASVs per bioclimatic region. This means that these ASVs are found in only one bioclimatic region 
prop_unique_MR<- endemic_MR_prev%>%
  group_by(MR)%>% 
  summarize(prop=prev/endemic_oneSPD_n)
colnames(prop_unique_MR) <- c("Koppen_clim","prop_unique")

## Filter the table of prevalence >1 & <10 SPDs (5% of sites)
twoandnine<- howmanyasv %>% filter(Count_nb > 1 & Count_nb < 10) ## 4717 ASVs 4717/15218= 31% of ASVs found in 1 to 9 SPDs
twoandnine_n <- twoandnine %>% summarize(n())
twoandnine_n/length(unique(row.names(prune_asv16s_df)))

# merging twoandnine with ASV table to get sample id and for each sample we attribute a MR
merge_twoandnine <- merge(twoandnine, asvdfmelt, by.x="ASV",by.y="Var2")
merge_twoandnine$MR <- vapply((merge_twoandnine$Var1), function(x) metadata_nomis$koppen_clim[metadata_nomis$Sample == x], FUN.VALUE = character(1))

colnames(merge_twoandnine) <- c("ASV","prev_GFS","glname","nb_count","MR")

# How many times an ASV is identified in a mountain range, sum by gl_code
# we count the nb of lines by ASV and gl_code. 
twoandnine_end<- merge_twoandnine %>% group_by(MR,ASV) %>% summarize(prev=n()) 
twoandnine_end$Color <- "B_twoandnine"

# now we would like it to be unique in this MR! we mutate to count the nb of times this ASVs is observed (in how many MR),then we filter ASVs found in one MR
# then we group by MR and we sum by N to know how many ASVs per MR
endemism_twoandnine <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) 

endemism_twoandnine_plot <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

# Number of ASVs that are specific to a mountain range (present in 2-9 GFSs)
sum_endemic_twoandnine <- endemism_twoandnine %>% summarize(sum=sum(number))
sum_endemic_twoandnine/length(unique(row.names(prune_asv16s_df)))

## Proportion by MR
prop_unique_twonine<- endemism_twoandnine%>%
  group_by(MR)%>% 
  summarize(prop=number/twoandnine_n)

colnames(prop_unique_twonine)<-c("Koppen_clim", "prop_unique")

tenandmore<- howmanyasv %>% filter(Count_nb >= 10) ## 1197 -> 7.9% of all ASVs
tenandmore_n <- tenandmore %>% summarize(n())
tenandmore_n/length(unique(row.names(prune_asv16s_df)))

# merge tenandmore with asv table to get the sample id and then for each sample we attribute a MR 
merge_tenandmore <- merge(tenandmore, asvdfmelt, by.x="ASV",by.y="Var2")
merge_tenandmore$MR <- vapply((merge_tenandmore$Var1), function(x) metadata_nomis$koppen_clim[metadata_nomis$Sample == x], FUN.VALUE = character(1))

colnames(merge_tenandmore) <- c("ASV","prev_GFS","glname","nb_count","MR")

# sum ASvs
tenandmore_end<- merge_tenandmore %>% group_by(MR,ASV) %>% summarize(prev=sum(nb_count>0))

tenandmore_end$Color <-"A_tenandemore"

# We would like it to be unique from this MR! mutate to count the nb of times when ASV has been identified, then filter only ASVs that are observed 1X in the MR
# group by MR and sum by N to know how many ASVs per MR

endemism_tenandmore <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 329/38983=0.8% (>=10 GFSs).
###2024 -- 648/54019=1.2% of total ASVs are specific to this section (>=10GFs)

endemism_tenandmore_plot <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

sum(endemism_tenandmore$number)

## Control Sanity check! ##
control_asv <- asvdfmelt
control_asv$MR <- vapply((control_asv$Var1), function(x) metadata_nomis$koppen_clim[metadata_nomis$Sample == x], FUN.VALUE = character(1))
## We want to know how many ASVs are specific, and found in only one bioclimatic region (hereafter MR)! We start from the original ASV table

controlasv_end<- control_asv %>% group_by(MR,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n())

# A tibble: 3 × 2
#MR                           number
#<chr>                        <int>
# 1 Continental cold/Tundra   2769
# 2 Continental hot summer    3252
# 3 Temperate                 4869

#10890 that's the sum. 10890/15218=71.6%

#write.csv(controlasv_end,"202403_endemic_list.csv")

## Number of Observed ASVs per bioclimatic region
controlasv_total<- control_asv %>% group_by(MR,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(MR)%>%summarize(sumii=sum(n()))


##Data 2024
# A tibble: 3 × 2
#MR                           number  sumi i
#<chr>                        <int>
# 1 Continental cold/Tundra   2769    5920
# 2 Continental hot summer    3252    6943
# 3 Temperate                 4869    8525

## Proportion of specific ASVs per mountain range
prop_endemic_MR <- merge(controlasv_end, controlasv_total, by="MR")
prop_endemic_MR$prop_endemic <- prop_endemic_MR$number/prop_endemic_MR$sumii

#MR                           number  sumii   prop_endemic
#1  Continental cold/Tundra   2769    5920    0.4677365
#2  Continental hot summer    3252    6943    0.4683854
#3               Temperate    4869    8525    0.5711437

## Proportion of unique ASVs -- 50.2% of specific ASVs are unique to one SPD!
endemic_oneSPD_n/sum(controlasv_end$number)

#n
#1 0.8543618

## Graphs
# We rbind the 3 different dataframes to plot the ASVs that are specific to the different SPD bioclimatic areas
endemic_oneSPD_plot$n <- 1
df_full <- rbind(endemic_oneSPD_plot, endemism_tenandmore_plot, endemism_twoandnine_plot)
niveaux <- c("Temperate","Continental cold/Tundra","Continental hot summer")
df_full$MR<- factor(df_full$MR, levels = niveaux)
df_full$MR <- fct_rev(df_full$MR)

ggplot(df_full, aes(fill=Color, y=MR)) + 
  geom_bar(position="stack", stat="count")+theme_minimal()

#proportion of ASVs that are endemic and specific to SPDs per bioclimatic area
prop_end_spe <- df_full %>% group_by(MR, Color) %>% summarize(sumii=sum(n()))
prop_end_spe2 <- df_full %>% group_by(MR) %>% summarize(sumii=sum(n()))

## Number of Observed ASVs per bioclimatic area
asv_total<- control_asv %>% group_by(MR,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(MR)%>%summarize(sumii=sum(n()))
asv_total <- inner_join(prop_end_spe2, asv_total, by = 'MR')
asv_total$sumii <- asv_total$sumii.y-asv_total$sumii.x
asv_total$Color <- "A_Other"
asv_total <- asv_total[,c(1,5,4)]


## Merge data
all_prop <- rbind(prop_end_spe, asv_total)

ggplot(all_prop, aes(x = MR, y = sumii, fill = Color)) +
  geom_bar(stat = "identity") +
  labs(x = "Climate Type", y = "Sum", title = "Stacked Bar Plot by Climate Type and Color") +
  theme_minimal() + coord_flip() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Summarize total and relevant colors per MR
prop_df <- all_prop %>%
  group_by(MR) %>%
  summarise(
    total = sum(sumii),
    highlight_sum = sum(sumii[Color %in% c("B_twoandnine", "C_uniquetoone")]),
    proportion = highlight_sum / total
  )

# Main stacked bar plot
ggplot(all_prop, aes(x = MR, y = sumii, fill = Color)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(data = prop_df, 
            aes(x = MR, y = total + 500,  # Offset above the tallest bar
                label = paste0(round(proportion * 100, 1), "%")),
            inherit.aes = FALSE) +
  scale_y_continuous(
    breaks = seq(0, max(prop_df$total) + 500, by = 3000),  # Y-axis breaks every 1000
    expand = expansion(mult = c(0, 0.05))  # Avoid too much space below bars
  ) +
  labs(x = "Bioclimatic area", y = "Number of specific ASVs") +
  theme_bw() + coord_flip()+
  scale_fill_manual(values = c(
    "A_Other" = "gray70",
    "B_twoandnine" = "darkorange",
    "C_uniquetoone" = "steelblue"
  ))


