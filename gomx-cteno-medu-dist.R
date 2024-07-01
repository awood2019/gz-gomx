#making figures for data exploration & analysis of ctenophore and medusozoan paper

#------------make phyloseq object of GoMx ASVs with taxonomy for ctenophores and medusozoans----------

#load packages
library(dplyr)
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(dendextend) ; packageVersion("dendextend") 
library(DESeq2)
library(viridis) ; packageVersion("viridis") 
library(ggplot2)

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-cteno-medu-dist")

#read in ctenophore and medusozoan classifier results w/ taxonomy and ASV
cteno_tax <- read_delim("Cteno_classifier_results.tsv", show_col_types = FALSE)
medu_tax <- read_delim("Medu_classifier_results.tsv", show_col_types = FALSE)

#read in original phyloseq taxonomy table w/ all GoM ASVs, sequence, & classified taxonomy
all_tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) 

#rename seq column in classifier results
names(cteno_tax)[names(cteno_tax) == "seq"] <- "Sequence"
names(medu_tax)[names(medu_tax) == "seq"] <- "Sequence"

#join ctenophore and medusozoan classifier results by sequence
tax_cteno_medu <- full_join(cteno_tax, medu_tax)

#join ctenophore classifier results with taxonomy table from phyloseq object to make new taxonomy table of more detailed cteno/medu taxonomy
tax_tab_cteno_medu <- left_join(all_tax_tab, tax_cteno_medu, by = c("Sequence", "Phylum"))

#delete unnecessary columns
tax_tab_cteno_medu <- select(tax_tab_cteno_medu, -`unique abundance`)

#save table
write_csv(tax_tab_cteno_medu, "~/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_tab_cteno_medu.csv")

#load components of phyloseq object: taxonomy table, count table, and sample data table
tax_tab <- read_csv("tax_tab_cteno_medu.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & classified taxonomy
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_cls: includes taxonomic information for each representative (ASV) sequence
class <- tax_tab$Class #pulling out class column from taxonomy table
tax_tab_cls <- tibble(class) #making phyla into a tibble containing class for each sequence
tax_tab_cls <- as.matrix(tax_tab_cls) #make tibble into matrix
row.names(tax_tab_cls) <- tax_tab$Sequence #make sequence column the row names
tax_tab_cls[is.na(tax_tab_cls)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_cls: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_cls)
count_tab_cls <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_cls) <- count_tab$...1 #make sequences the row names (ignore warning message)

#sample_info_tab_cls: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab$sample_or_control <- ifelse(sample_info_tab$`CTD/ROV` == "Negative" | sample_info_tab$`CTD/ROV` == "NTC", "control", "sample")
sample_info_tab_cls <- sample_info_tab
sample_info_tab_cls <- sample_info_tab_cls[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_cls) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object 
ASV_physeq <- phyloseq(otu_table(count_tab_cls, taxa_are_rows = TRUE), tax_table(tax_tab_cls), sample_data)
ASV_physeq <- prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'cteno_medu_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe 
df_ASV_physeq <- ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe

#save phyloseq dataframe
write.csv(df_ASV_physeq, "~/Desktop/AW RStudio/results/ASV_physeq.csv")

#make phyloseq object with only ctenophores and medusozoans, excluding those that weren't assigned taxonomy
cm_physeq <- subset_taxa(ASV_physeq, class == "Alphaclass" |class == "Betaclass"|class == "Hydrozoa" |class == "Scyphozoa") 
df_cm_physeq <- cm_physeq %>% psmelt() 

#testing by comparing to old classifiers-------------
#read in ctenophore and medusozoan classifier results w/ taxonomy and ASV
old_cteno_tax <- read_delim("Cteno_classifier_results_old.tsv", show_col_types = FALSE)
old_medu_tax <- read_delim("Medu_classifier_results_old.tsv", show_col_types = FALSE)

#read in original phyloseq taxonomy table w/ all GoM ASVs, sequence, & classified taxonomy
all_tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) 

#rename seq column in classifier results
names(old_cteno_tax)[names(old_cteno_tax) == "seq"] <- "Sequence"
names(old_medu_tax)[names(old_medu_tax) == "seq"] <- "Sequence"

#join ctenophore and medusozoan classifier results by sequence
tax_cteno_medu <- full_join(old_cteno_tax, old_medu_tax)

#join ctenophore classifier results with taxonomy table from phyloseq object to make new taxonomy table of more detailed cteno/medu taxonomy
tax_tab_cteno_medu <- left_join(all_tax_tab, tax_cteno_medu, by = c("Sequence", "Phylum"))

#delete unnecessary columns
tax_tab_cteno_medu <- select(tax_tab_cteno_medu, -`unique abundance`)

#save table
write_csv(tax_tab_cteno_medu, "~/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_tab_cteno_medu.csv")

#load components of phyloseq object: taxonomy table, count table, and sample data table
tax_tab <- read_csv("tax_tab_cteno_medu.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & classified taxonomy
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_cls: includes taxonomic information for each representative (ASV) sequence
class <- tax_tab$Class #pulling out class column from taxonomy table
tax_tab_cls <- tibble(class) #making phyla into a tibble containing class for each sequence
tax_tab_cls <- as.matrix(tax_tab_cls) #make tibble into matrix
row.names(tax_tab_cls) <- tax_tab$Sequence #make sequence column the row names
tax_tab_cls[is.na(tax_tab_cls)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_cls: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_cls)
count_tab_cls <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_cls) <- count_tab$...1 #make sequences the row names (ignore warning message)

#sample_info_tab_cls: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab$sample_or_control <- ifelse(sample_info_tab$`CTD/ROV` == "Negative" | sample_info_tab$`CTD/ROV` == "NTC", "control", "sample")
sample_info_tab_cls <- sample_info_tab
sample_info_tab_cls <- sample_info_tab_cls[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_cls) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object 
old_ASV_physeq <- phyloseq(otu_table(count_tab_cls, taxa_are_rows = TRUE), tax_table(tax_tab_cls), sample_data)
old_ASV_physeq <- prune_taxa(taxa_sums(old_ASV_physeq) > 0, old_ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'cteno_medu_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe 
olddf_ASV_physeq <- old_ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe

df_ASV_physeq <- dplyr::rename(df_ASV_physeq, c(class2 = "class"))
test <- left_join(olddf_ASV_physeq, df_ASV_physeq)

#-------------------------decontam to remove ASVs in blanks--------------------

#load packages needed for decontam
library(decontam)

#how many unique ASVs are in cm_physeq to start?
length(unique(df_cm_physeq$OTU)) #459

#inspect library sizes (# of reads per sample) to see if sample was a true pos sample or neg control
df <- as.data.frame(sample_data(cm_physeq)) #put sample data into data frame for ggplot
df$LibrarySize <- sample_sums(cm_physeq) #library size = total # of ASVs in each sample
df$Index <- seq(nrow(df)) #index = row of df
ggplot(data = df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point() #graph of reads per sample

#identify contaminants
sample_data(cm_physeq)$is.neg <- sample_data(cm_physeq)$sample_or_control == "control"
contamdf.prev <- isContaminant(cm_physeq, method="prevalence", neg="is.neg") #identifies contaminants by comparing prevalence of each ASV in samples to prevalence in controls 
table(contamdf.prev$contaminant) #0 contaminants, 459 true samples

#more agressive contaminant threshhold
contamdf.prev05 <- isContaminant(cm_physeq, method = "prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant) #4 contaminants, 473 true samples

#look at # of times these taxa were observed in neg controls and pos samples
#phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(cm_physeq, function(abund) 1*(abund>0)) 
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "sample", ps.pa)

#dataframe of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#---------------UniFrac distance matrix of samples for NMDS plot---------------

#load packages needed for NMDS ordination plot
library(vegan)
library(ape)
library(Biostrings)
library(tidysq)
library(ggalt)
library(phytools)

#take out all negative and blank samples
cm_physeq <- subset_samples(cm_physeq, CTD.ROV != "Negative" & CTD.ROV !="NTC")

#remove samples with 0 organisms so the distance matrix can be calculated
cm_physeq <- prune_samples(sample_sums(cm_physeq) >0, cm_physeq) 
df_cm_physeq <- cm_physeq %>% psmelt()

#combine with taxonomy table to count taxa
total_taxa <- left_join(df_cm_physeq, tax_cteno_medu, by = c("OTU" = "Sequence"))
length(unique(total_taxa$Family)) #35 taxa total

#export table of ASVs so we can put it into Geneious to make a phylogenetic tree
tax_tab_cm <- tax_table(cm_physeq) #take component of ASVs out of phyloseq object
write.csv(tax_tab_cm, "~/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_cteno_medu.csv") #save as csv file

#add tree of ctenophores and medusozoans to phyloseq object
cm_tree <- read.nexus("/Users/quattrinia/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_cteno_medu_updated_clsfr alignment FastTree Tree.nex") #import Geneious tree of aligned ctenos and medus
rooted_tree <- midpoint.root(cm_tree)
cm_physeq_tree <- merge_phyloseq(cm_physeq, rooted_tree) #merge tree with existing phyloseq object

#calculate Unifrac distance matrix 
uni_matrix <- UniFrac(cm_physeq_tree, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #unweighted UniFrac
uni_matrix_wt <- UniFrac(cm_physeq_tree, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE) #weighted UniFrac

#--------------adding data for temperature, dissolved oxygen, and pH to sample metadata table for dbRDA--------------

#load packages
library(stringr)

#import chemistry metadata table
chem_metadata <- read_csv("PS21_eDNA_sample_data_updated.csv")

#create new column in sample metadata table for SH number
sample_info_tab[c('SH number', 'meso', 'octo')] <- str_split_fixed(sample_info_tab$`Sample name`, "-", 3)
sample_info_tab <- select(sample_info_tab, -`meso` & -`octo`) #delete unnecessary columns

#join chemistry metadata table and sample metadata table by SH number
metadata <- left_join(sample_info_tab, chem_metadata, by = "SH number")

#take out all irrelevant columns added from chemistry metadata table
metadata <- select(metadata, -`Cast/Dive #` & -`Date` & -`Site.y` & -`Niskin.y` & -`Cruise sample ID` & -`Latitude` & -`Longitude` & -`Time fired` & -`Target depth` & -`Actual trip depth` & -`Altitude.y` & -`Volume.y` & -`Nutrients` & -`DOC` & -`Cell counts`)

#rename columns to not have spaces
metadata <- dplyr::rename(metadata, c(CTD.Salinity = "CTD Salinity", ROV.Salinity = "ROV Salinity", CTD.Temp = "CTD Temperature", ROV.Temp = "ROV Temperature", CTD.pH = "CTD pH", Lab.pH = "Lab pH"))
metadata <- dplyr::rename(metadata, c(Site = "Site.x", Niskin = "Niskin.x", Altitude = "Altitude.x", Volume = "Volume.x"))

#take out blanks
metadata <- metadata %>%
  subset(!(`CTD/ROV` == "Negative" | `CTD/ROV` == "NTC"))

#save metadata table so we can combine separate temp, salinity, pH variable columns for CTD/ROV
write.csv(metadata, "dbrda_metadata.csv")

#import updated metadata table with combined variable columns for CTD and ROV
metadata_updated <- read_csv("dbrda_metadata_updated.csv", show_col_types = FALSE)

#code CTD/ROV columns as dummy variables
metadata_updated$CTD <- as.numeric(ifelse(metadata_updated$`CTD/ROV` == "CTD", "1", "0"))
metadata_updated$ROV <- as.numeric(ifelse(metadata_updated$`CTD/ROV` == "ROV", "1", "0"))

#-------------normalize environmental variables for dbRDA--------------

#copy metadata table containing variables we're going to normalize
norm_metadata <- select(metadata_updated, c(`Site`, `Sample`, `Sample name`, `SH number`, `CTD/ROV`, `Depth`, `Altitude`, `Salinity`, `Temp`, `pH`, `CTD`, `ROV`)) 

#normalize environmental variables
norm_metadata$Depth <- norm_metadata %>% select(Depth) %>% decostand("range", na.rm = TRUE)
norm_metadata$Altitude <- norm_metadata %>% select(Altitude) %>% decostand("range", na.rm = TRUE)
norm_metadata$Salinity <- norm_metadata %>% select(Salinity) %>% decostand("range", na.rm = TRUE)
norm_metadata$Temp <- norm_metadata %>% select(Temp) %>% decostand("range", na.rm = TRUE)
norm_metadata$pH <- norm_metadata %>% select(pH) %>% decostand("range", na.rm = TRUE)

norm_metadata$Depth <- unlist(norm_metadata$Depth)
norm_metadata$Altitude <- unlist(norm_metadata$Altitude)
norm_metadata$Salinity <- unlist(norm_metadata$Salinity)
norm_metadata$Temp <- unlist(norm_metadata$Temp)
norm_metadata$pH <- unlist(norm_metadata$pH)

#-----------------------dbRDA and PERMANOVA for environmental variables---------------------

#load packages
library(ggordiplots)
library(BiodiversityR)
library(ggrepel)
library(ggpubr)

#perform dbRDA on unweighted UniFrac distance matrix
dbrda <- dbrda(formula = uni_matrix ~ Depth + Altitude + Salinity + Temp + pH + ROV + CTD, data = norm_metadata) 

#view unweighted dbRDA results
print(dbrda)
dbrda_sum <- summary(dbrda)
dbrda_sum$concont #constrained eigenvalues show proportion of the variance explained by these explanatory variables
dbrda_sum$biplot #coefficients for each explanatory variable on each axis
0.3807 * 0.3855 #dbRDA1 explains 38% of the overall 39% of variation accounted for by all variables, so 14.7% of the total variation
0.2966 * 0.3855 #dbRDA2 explains 11.4% of total variation
ordiplot <- ordiplot(dbrda) #quick visualization

#run PERMANOVA on unweighted dbRDA
pnova <- adonis2(formula = uni_matrix ~ Depth + Altitude + Salinity + Temp + pH + ROV + CTD, data = norm_metadata)
print(pnova) #shows variables that significantly affect the variation and what % of the variation they account for

#perform dbRDA on weighted UniFrac distance matrix
dbrda_wt <- dbrda(formula = uni_matrix_wt ~ Depth + Altitude + Salinity + Temp + pH + ROV + CTD, data = norm_metadata) #weighted UniFrac

#view weighted dbRDA results
print(dbrda_wt)
dbrda_sum_wt <- summary(dbrda_wt)
dbrda_sum_wt$concont #constrained eigenvalues show proportion of the variance explained by these explanatory variables
dbrda_sum_wt$biplot #coefficients for each explanatory variable on each axis
0.6455 * 0.4578 #dbRDA1 explains 65% of the overall 45.8% of variation accounted for by all variables, so 30% of the total variation
0.1634 * 0.4578 #dbRDA2 explains 7.5% of overall variation
ordiplot_wt <- ordiplot(dbrda_wt) #quick visualization

#run PERMANOVA on weighted dbRDA
pnova_wt <- adonis2(formula=uni_matrix_wt ~ Depth + Altitude + Salinity + Temp + pH + ROV + CTD, data = norm_metadata)
print(pnova_wt)

#plot unweighted dbRDA
#make dbrda and ordiplot info ggplot friendly so we can plot it
sites.long <- sites.long(ordiplot, env.data=norm_metadata) #extract location info of sites in ordiplot & add environmental data
sites.long
axis.long <- axis.long(dbrda, choices=c(1,2)) #extract axis label info 
axis.long
arrows.scale <- ordiplot$biplot #extract arrow coordinate info from dbrda plot
CTD <- as.numeric(list('-0.73678729','0.32738800'))
arrows.scale <- rbind(arrows.scale, CTD)
arrows.scale <- as.data.frame(arrows.scale)
arrows.scale$label <- rownames(arrows.scale)
unwt_plot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long, 
             aes(x=axis1, y=axis2, colour=Site, shape=CTD.ROV), 
             size=3) +
  geom_segment(data=arrows.scale,
               aes(x=0, y=0, xend=dbRDA1, yend=dbRDA2),
               colour="black", size=0.5, arrow=arrow(length=unit(0.1, "inches"))) +
  geom_text_repel(data=arrows.scale,
                  aes(x=dbRDA1*1.3, y=dbRDA2*1.3, label=label),
                  size=4) +
  labs(shape = 'Sampling Method') +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1) +
  theme_classic()
unwt_plot

#plot weighted dbRDA
#make dbrda and ordiplot info ggplot friendly so we can plot it
sites.long.wt <- sites.long(ordiplot_wt, env.data=norm_metadata) #extract location info of sites in ordiplot & add environmental data
axis.long.wt <- axis.long(dbrda_wt, choices=c(1,2)) #extract axis label info 
arrows.scale.wt <- ordiplot_wt$biplot #extract arrow coordinate info from dbrda plot
CTD <- as.numeric(list('0.9205298', '-0.1159485'))
arrows.scale.wt <- rbind(arrows.scale.wt, CTD)
arrows.scale.wt <- as.data.frame(arrows.scale.wt)
arrows.scale.wt$label <- rownames(arrows.scale.wt)
wt_plot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long.wt[1, "label"]) +
  ylab(axis.long.wt[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long.wt, 
             aes(x=axis1, y=axis2, colour=Site, shape=CTD.ROV), 
             size=3) +
  geom_segment(data=arrows.scale.wt,
               aes(x=0, y=0, xend=dbRDA1, yend=dbRDA2),
               colour="black", size=0.5, arrow=arrow(length=unit(0.1, "inches"))) +
  geom_text_repel(data=arrows.scale.wt,
                  aes(x=dbRDA1*1.3, y=dbRDA2*1.3, label=label),
                  size=4) +
  labs(shape='Sampling Method') +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1) +
  theme_classic()
wt_plot

#plot unweighted and weighted dbRDA together
dbrda_plot <- ggarrange(unwt_plot, wt_plot,
                     ncol = 2, nrow = 1)
dbrda_plot


#trying ggplot with envfit?
envfit <- envfit(ordiplot, env=norm_metadata) #envfit calculates significant environmental variables
envfit
data.envfit <- data.frame(r=envfit$vectors$r, p=envfit$vectors$pvals)
vectors.envfit <- data.frame(envfit$vectors$arrows)
vectors.long <- vectorfit.long(envfit)
vectors.long
vectors.scale <- vectors.long[c(1:7),c(2:3)] * vectors.long$r
ggplot() + arrows(0,0, vectors.scale[,1], vectors.scale[,2], len=0.1)
ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long, 
             aes(x=axis1, y=axis2, colour=Site, shape=CTD.ROV), 
             size=3) +
  geom_segment(data=vectors.scale,
               aes(x=0, y=0, xend=axis1, yend=axis2),
               colour="black", size=0.5, arrow=arrow(length=unit(0.1, "inches"))) +
  geom_text_repel(data=vectors.scale,
                  aes(x=axis1*1.3, y=axis2*1.3, label=label),
                  size=4)+
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1) +
  theme_classic()

#plot dbRDA using gg_envfit, plots ordination with env variables as vectors colored by site
envfit <- gg_envfit(ord = dbrda, env = norm_metadata, groups = norm_metadata$Site, arrow.col = "black", len = 0.2, pt.size = 2) #displays env variables as vectors
names(envfit) #list components of envfit object
envfit.data <- envfit$df_ord #capture the df of point coordinates for ordination
envfit.data$Sampling_Method <- norm_metadata$`CTD/ROV`#add sampling method to ordination df
envfit$df_ord <- envfit.data #put modified df of point coordinates back into envfit object
envfit.plot <- envfit$plot #capture the plot
envfit.plot.df <- envfit.plot$data #capture the df of point coordinates for ordination from plot object
envfit.plot.df$Sampling_Method <- norm_metadata$`CTD/ROV` #add sampling method column to envfit plot df
envfit.plot$data <- envfit.plot.df #replace original df of point coordinates w/ updated one
envfit.plot + theme_classic() + labs(color = "Site", shape = "Sampling Method", x = "dbRDA1 (17.9%)", y = "dbRDA2 (12%)") +
  geom_point(data = envfit.plot.df, aes(x=x, y=y, shape = Sampling_Method), size = 3) 
vegan::ord

#plot dbRDA using gg_envfit, plots ordination with env variables as vectors colored by site
envfit_wt <- gg_envfit(ord = dbrda_wt, env = norm_metadata, groups = norm_metadata$Site, arrow.col = "black", len = 0.2, pt.size = 2) #displays env variables as vectors
names(envfit_wt) #list components of envfit object
envfit.data.wt <- envfit_wt$df_ord #capture the df of point coordinates for ordination
envfit.data.wt$Sampling_Method <- norm_metadata$`CTD/ROV`#add sampling method to ordination df
envfit_wt$df_ord <- envfit.data.wt #put modified df of point coordinates back into envfit object
envfit.plot.wt <- envfit_wt$plot #capture the plot
envfit.plot.df.wt <- envfit.plot.wt$data #capture the df of point coordinates for ordination from plot object
envfit.plot.df.wt$Sampling_Method <- norm_metadata$`CTD/ROV` #add sampling method column to envfit plot df
envfit.plot.wt$data <- envfit.plot.df.wt #replace original df of point coordinates w/ updated one
envfit.plot.wt + theme_bw() + labs(shape = "Sampling Method", x = "dbRDA1 (21.0%)", y = "dbRDA2 (8.0%)") +
  geom_point(data = envfit.plot.df.wt, aes(x=x, y=y, shape = Sampling_Method), size = 3) 

#---------------------------alpha diversity across sites for ROV------------------

#load packages
library(picante)
library(genefilter)
library(btools)
library(ggpubr)

#prune ASVs classified as species that aren't observed 
alphdiv_physeq <- prune_taxa(taxa_sums(cm_physeq) > 0, cm_physeq)

#estimate alpha diversity measures
richness <- estimate_richness(cm_physeq)
head(richness)

#specify site order by increasing depth in vector so it can be added to ggplot
site_order <- c('Stetson Bank', 'Bright Bank', 'Viosca Knoll', 'Green Canyon')

#plot taxa richness and Shannon diversity
richness <- plot_richness(cm_physeq, x="Site", color = "Site", measures = c("Chao1", "Shannon")) +
  geom_boxplot(alpha=0.6, aes(x=factor(Site, levels = site_order))) +
  ggsci::scale_colour_npg() +
  theme_classic()
richness

#calculate Faith's phylogenetic diversity
pd.result <- estimate_pd(cm_physeq_tree)

#merge Faith's phylogenetic diversity results w/ sample data to get site info
sample_data_cm <- cm_physeq_tree@sam_data #take out sample data table from physeq object
pd.result$File.name <- row.names(pd.result)
pd_sampledata <- left_join(pd.result, sample_data_cm)

#boxplot of Faith's phylogenetic diversity
boxplot(PD~Site, data=pd_sampledata)
pd <- ggplot(pd_sampledata, aes(x=Site, y=PD, color = Site)) +
  ylab("Phylogenetic Diversity Measure") +
  geom_boxplot(alpha=0.6, aes(x=factor(Site, levels=site_order))) +
  geom_point(alpha=0.6) +
  ggsci::scale_colour_npg() +
  theme_classic()
pd

#put Chao1, Shannon, and phylogenetic diversity in same plot
alphdiv <- ggarrange(richness, pd,
                     ncol = 1, nrow = 2)
alphdiv

#--------------test if we can combine bottom CTDs with ROVs for benthic vs pelagic comparison-------------

#filter out only ROV and bottom CTD samples for physeq object
bottom_samples <- cm_physeq %>%
  subset_samples(Site == "Bright Bank" | Site == "Viosca Knoll")
bottom_samples <- bottom_samples %>%
  subset_samples(Altitude < 5.8)

#filter out only ROV and bottom CTD samples for metadata table
bottom_metadata <- metadata %>%
  subset(Site == "Bright Bank" | Site == "Viosca Knoll")
bottom_metadata <- bottom_metadata %>%
  subset(Altitude < 5.8)

#export table of ASVs so we can put it into Geneious to make a phylogenetic tree
tax_bottom <- tax_table(bottom_samples) #take component of ASVs out of phyloseq object
write.csv(tax_bottom, "~/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_bottom.csv") #save as csv file

#add tree of ctenophores to each phyloseq object
bottom_tree <- read.nexus("/Users/quattrinia/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_bottom alignment FastTree Tree.nex") #import Geneious tree of aligned ctenophores
bottom_rooted_tree <- midpoint.root(bottom_tree)
bottom_samples_tree <- merge_phyloseq(bottom_samples, bottom_rooted_tree) #merge tree with existing phyloseq object for abundance >0

#compute distance matrix for bottom_samples
bottom_matrix <- UniFrac(bottom_samples_tree, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) 

#PERMANOVA for bottom samples
pnova_bottom <- adonis2(formula = bottom_matrix ~ `CTD/ROV`, data = bottom_metadata)
print(pnova_bottom)

#-------------------plot for benthic comparison across all sites-----------------

#since all ROVs are benthic, can plot by CTD/ROV variable to represent benthic samples

#load packages
library(gridExtra)

#create table from metadata table with only variables of interest
benthic_table <- cm_physeq %>% psmelt()
benthic_table <- benthic_table %>%
  select(OTU, Abundance, Site, Depth, Altitude, CTD.ROV, depth_bin, class)

#select only ASVs collected at the seafloor 
benthic_table <- benthic_table %>% subset(CTD.ROV == "ROV")

#select only ASVs with abundance > 0
benthic_table <- benthic_table %>% filter(Abundance > 0)

#join metadata table with taxonomy table
benthic_tax_table <- left_join(benthic_table, tax_cteno_medu, by =c("OTU" = "Sequence"))

#take all ASVs not classified to family level out
fam_benthic <- benthic_tax_table %>% subset(!is.na(Family))
write.csv(fam_benthic, "~/Desktop/AW RStudio/results/gomx-cteno-medu-dist/benthic_dist.csv") #save table as csv so we can analyze in Excel

#plot benthic communities with separate legends 
grid_labs <- c("Medusozoa", "Ctenophora") #make vector for modifying grid names
names(grid_labs) <- c("Cnidaria", "Ctenophora") #add original grid names to vector
separate.plots <- by(data=fam_benthic, INDICES = fam_benthic$Phylum, FUN = function(m) {
  m <- droplevels(m)
  m <- ggplot(m, aes(x=factor(Site, level = c("Stetson Bank", "Bright Bank", "Viosca Knoll", "Green Canyon")), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by family
    geom_bar(position="fill",stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
    facet_grid(.~Phylum, labeller = labeller(Phylum = grid_labs), scale = "free_x", space = "free_x") +
    scale_y_continuous(name = "Percentage of sequencing reads recovered",
                       labels = scales::percent) +
    xlab("Site") +
    theme_classic()+
    scale_fill_viridis(discrete=TRUE, option="plasma")
})
do.call(grid.arrange,separate.plots)

length(unique(fam_benthic$Family)) #number of unique families total - 23

#abundance of unique families at each site
sb <- fam_benthic %>% filter(Site == "Stetson Bank")
ggplot(sb, aes(x=factor(Family), y = Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90))

bb <- fam_benthic %>% filter(Site == "Bright Bank")
ggplot(bb, aes(x=factor(Family), y = Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90))

vk <- fam_benthic %>% filter(Site == "Viosca Knoll")
ggplot(vk, aes(x=factor(Family), y = Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90))

gc <- fam_benthic %>% filter(Site == "Green Canyon")
ggplot(gc, aes(x=factor(Family), y = Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90))

#---------------plot for pelagic comparison---------------

#create table from metadata table with only variables of interest
pelagic_table <- cm_physeq %>% psmelt()
pelagic_table <- pelagic_table %>%
  select(OTU, Abundance, Site, Depth, Altitude, CTD.ROV, depth_bin, class)

#select only ASVs collected at the water column 
pelagic_table <- pelagic_table %>% subset(CTD.ROV == "CTD")

#select only ASVs with abundance > 0
pelagic_table <- pelagic_table %>% filter(Abundance > 0)

#join metadata table with taxonomy table
pelagic_tax_table <- left_join(pelagic_table, tax_cteno_medu, by =c("OTU" = "Sequence"))

#take all ASVs not classified to family level out
fam_pelagic <- pelagic_tax_table %>% subset(!is.na(Family))
write.csv(fam_pelagic, "~/Desktop/AW RStudio/results/gomx-cteno-medu-dist/pelagic_dist.csv") #save table as csv so we can analyze in Excel

#plot pelagic communities with separate legends 
grid_labs <- c("Medusozoa", "Ctenophora") #make vector for modifying grid names
names(grid_labs) <- c("Cnidaria", "Ctenophora") #add original grid names to vector
separate.plots <- by(data=fam_pelagic, INDICES = fam_pelagic$Phylum, FUN = function(m) {
  m <- droplevels(m)
  m <- ggplot(m, aes(x=factor(Site, level = c("Stetson Bank", "Bright Bank", "Viosca Knoll", "Green Canyon")), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by family
    geom_bar(position="fill",stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
    facet_grid(.~Phylum, labeller = labeller(Phylum = grid_labs), scale = "free_x", space = "free_x") +
    scale_y_continuous(name = "Percentage of sequencing reads recovered",
                       labels = scales::percent) +
    xlab("Site") +
    theme_classic()+
    scale_fill_viridis(discrete=TRUE, option="plasma")
})
do.call(grid.arrange,separate.plots)

length(unique(fam_pelagic$Family)) #number of unique families total - 28

#-------------------plot for CTD vs ROV comparison------------------

#create table from metadata table with only variables of interest
depth_table <- cm_physeq %>% psmelt()
depth_table <- depth_table %>%
  select(OTU, Abundance, Site, Depth, Altitude, CTD.ROV, depth_bin, class)

#select only ASVs with abundance > 0
depth_table <- depth_table %>% filter(Abundance > 0)

#join metadata table with taxonomy table
depth_tax_table <- left_join(depth_table, tax_cteno_medu, by =c("OTU" = "Sequence"))

#figure out which depth ranges were sampled at so we can make depth bins for plots
unique(depth_tax_table$depth_bin)

#take all ASVs not classified to family level out
fam_depth <- depth_tax_table %>%
  subset(!is.na(Family))

#plot by family
ggplot(fam_depth, aes(x=factor(Depth), y = Abundance, fill = Family, group = factor(Site))) + #x-axis = depth, y-axis = ASV abundance - plotted by family
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(~CTD.ROV, scale = "free_x", space = "free_x") +
  scale_y_continuous(name = "Percentage of sequencing reads recovered",
                     labels = scales::percent) +
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="magma")

#abundance of ASVs by family
ggplot(fam_depth, aes(x=factor(Family), y=Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

#plots by sampling method
ctd <- fam_depth %>% filter(CTD.ROV == "CTD")
ggplot(ctd, aes(x=factor(Family), y=Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

rov <- fam_depth %>% filter(CTD.ROV == "ROV")
ggplot(rov, aes(x=factor(Family), y=Abundance)) +
  geom_col() +
  ylab("ASV Abundance") +
  xlab("Family") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
ggplot(rov, aes(x=factor(Depth), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by family
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#--------NOT DOING ANYMORE plot for mesophotic vs deep comparison--------

#subset depth_tax_table to only be depths below 40m
photic_tax_table <- depth_tax_table %>%
  subset(Depth > 40)

#change depth bins to 40-200 and >200
photic_tax_table$zone_bin <- ifelse(photic_tax_table$depth_bin == "(40,50]"|photic_tax_table$depth_bin == "(50,60]" |photic_tax_table$depth_bin == "(60,70]" |photic_tax_table$depth_bin == "(70,80]" |photic_tax_table$depth_bin == "(80,90]" |photic_tax_table$depth_bin == "(110,120]",
                                  "meso", "deep")

#take all the ASVs not classified to family level out of dataset
fam_photic <- photic_tax_table %>%
  subset(!is.na(Family))
                                   
#plot by family
ggplot(fam_photic, aes(x=factor(zone_bin, level=c('meso', 'deep')), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by family
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#take all the ASVs not classified to order level out of dataset
order_zone <- photic_tax_table %>%
  subset(!is.na(Order))

#plot by order
ggplot(order_zone, aes(x=factor(zone_bin, level=c('meso', 'deep')), y = Abundance, fill = Order)) + #x-axis = depth, y-axis = ASV abundance - plotted by family
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 




#-------------NOT DOING ANYMORE HIERACHICAL CLUSTERING, JUST KEEPING CODE HANDY----------------

#load packages
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)

#subset dataset to just ctenophores and medusozoans
df_cteno_medu <- df_ASV_physeq %>%
  subset(class == "Tentaculata" | class == "Nuda" | class == "Hydrozoa" | class == "Scyphozoa") %>%
  filter(Abundance >0) %>%
  filter(CTD.ROV != "Negative" & CTD.ROV != "NTC")
  
#see if there are any NA's in data and remove if there are
sum(is.na(df_cteno_medu)) #count total NA's
sum(is.na(uni_matrix))

#figure out which columns are numeric so we can scale them
as.data.frame(sapply(df_cteno_medu, class))

#scale the data so it's standardized and variables can be made comparable
df_cteno_medu$Abundance <- scale(df_cteno_medu$Abundance, df_cteno_medu$Marker)

#----agglomerative hierarchical clustering with hclust----
hc1 <- hclust(uni_matrix, method="complete") #perform clustering using UniFrac dissimilarity matrix
plot(hc1, cex = 0.6, hang = -1) #plot the dendrogram of clustered samples

#color-code clusters by cutting dendrogram into subgroups
sub_grps <- cutree(hc1, k=3) #cut tree into 3 groups
sub_grps #see which samples were put into which subgroup
col <- c("red2", "green4", "mediumblue") #define vector of colors 
col[sub_grps] #assign a color to each point in the cluster so clusters are color coded

#plot NMDS with hc1 results
nmdsplot_hc1 <- plot_ordination(cm_physeq_tree, nmds, shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  theme_classic()+
  geom_point(size=3, color = col[sub_grps])+
  scale_color_continuous(type="viridis", option="H", direction=-1) +
  ggtitle("NMDS of Ctenophore and Medusozoan ASVs using hclust")
nmdsplot_hc1

#visualize hc1 results in scatter plot
fviz_cluster(list(data = uni_matrix, cluster=sub_grps, labelsize = 2)) 

#----agglomerative hierarchical clustering with agnes----
hc2 <- agnes(uni_matrix, method="complete") #perform clustering
hc2$ac #find agglomerative coefficient - the higher it is, the stronger the clustering structure 
pltree(hc2, cex = 0.6, hang = -1) #plot dendrogram

#determining which clustering method is the best
m <- c("average", "single", "complete", "ward") #define vector of method names
names(m) <- c("average", "single", "complete", "ward") #set the names of the vector
ac <- function(x) {agnes(uni_matrix, method = x)$ac} #function to compute agglomerative coefficient
  
map_dbl(m, ac) #calculate coefficients using each clustering method

#using ward method because it had the highest coefficient
hc3 <- agnes(uni_matrix, method="ward") #perform clustering
pltree(hc3, cex=0.6, hang=-1) #plot dendrogram of clustering results

#color-code clusters by cutting dendrogram into subgroups
sub_grps <- cutree(hc3, k=4)#cut tree into 3 subgroups
table(sub_grps)#see how many points are in each subgroup
sub_grps <- cutree(hc3, k=3) 
table(sub_grps) 
sub_grps #see which samples were put into each subgroup
col <- c("red2", "green4", "mediumblue", "darkgoldenrod") #define vector of colors
col[sub_grps] #assign a color to each cluster so they're color coded

#plot NMDS plot
nmdsplot_hc3 <- plot_ordination(cm_physeq_tree, nmds, shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  theme_classic()+
  geom_point(size=3, color = col[sub_grps])+
  scale_color_continuous(type="viridis", option="H", direction=-1) +
  ggtitle("NMDS of Ctenophore and Medusozoan Samples using agnes and Ward method")
nmdsplot_hc3

#visualize hc3 results in scatterplot
fviz_cluster(list(data = uni_matrix, cluster=sub_grps, labelsize = 2)) 

#----divisive hierarchical clustering with diana----
hc4 <- diana(uni_matrix) #perform clustering
hc4$dc #find divisive coefficient
pltree(hc4, cex=0.6, hang = -1) #plot dendrogram

#color-code clusters by cutting dendrogram into subgroups
sub_grps <- cutree(hc3, k=3) #cut tree into 3 subgroups
table(sub_grps)
col <- c("red2", "green4", "mediumblue")
col[sub_grps]

#plot NMDS plot
nmdsplot_hc4<- plot_ordination(cm_physeq_tree, nmds, shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  theme_classic()+
  geom_point(size=3, color = col[sub_grps])+
  scale_color_continuous(type="viridis", option="H", direction=-1) +
  ggtitle("NMDS of Ctenophore and Medusozoan Samples using diana")
nmdsplot_hc4

#visualize hc4 results in scatterplot
fviz_cluster(list(data = uni_matrix, cluster=sub_grps, labelsize = 2)) 

#----NOT DOING ANYMORE PCOA, JUST KEEPING THE CODE HANDY------ 
#perform PCoA ordination, taking into account relative abundance
cteno_medu_pcoa <- ordinate(cm_physeq_tree, method = "PCoA", uni_matrix, weighted=TRUE) #??idk why but Cailliez correction isn't working

#make PCOA plot for phyloseq object
plot_ordination(cm_physeq_tree, cteno_medu_pcoa, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Ctenophore and Medusozoan ASVs with Abundance >0")

#plotting distance matrix vs points calculated by PCoA to see how well ordination fits the data
plot(uni_matrix, dist(cteno_medu_pcoa$vectors), xlab = "UniFrac distance", ylab = "Distance mapped by PCoA") 

#calculating % of variance explained by PCoA axes
pcoa_var <- round(cteno_medu_pcoa$values$Eigenvalues*100/sum(cteno_medu_pcoa$values$Eigenvalues),1)

#-----NOT DOING ANYMORE Unifrac of unmerged sample replicates--------
#take out all negative and blank samples
cm_physeq <- subset_samples(cm_physeq, CTD.ROV != "Negative" & CTD.ROV !="NTC")
df_cm_physeq <- cm_physeq %>% psmelt()

#remove rows with 0 organisms so the distance matrix can be calculated
cm_physeq <- prune_samples(sample_sums(cm_physeq) >0, cm_physeq) 

#add tree of ctenophores and medusozoans to phyloseq object
cm_tree <- read.nexus("/Users/quattrinia/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_cteno_medu alignment FastTree Tree.nex") #import Geneious tree of aligned ctenos and medus
rooted_tree <- midpoint.root(cm_tree)
cm_physeq_tree <- merge_phyloseq(cm_physeq, rooted_tree) #merge tree with existing phyloseq object

#-------NOT DOING ANYMORE nmds of phyloseq object---------
#perform NMDS ordination for phyloseq object
nmds <- ordinate(cm_physeq_tree, method = "NMDS", uni_matrix, weighted=TRUE) 

#make NMDS plot for phyloseq object
plot_ordination(cm_physeq_tree, nmds, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  theme_classic()+
  geom_point(size=3)+
  scale_color_continuous(type="viridis", option="H", direction=-1) +
  ggtitle("NMDS of Ctenophore and Medusozoan Samples")

#plotting distance matrix vs points calculated by NMDS
plot(uni_matrix, dist(nmds$points), xlab = "UniFrac distance", ylab = "Distance mapped by NMDS") 

#to check, Shepard plots to determine how well NMDS fit the data
stressplot(nmds, uni_matrix)

#------NOT DOING ANYMORE relative abundance calculations--------
#unifrac calcs relative abundance for you so don't need this

#convert sample data table columns to factor so we can merge samples
df <- as.data.frame(lapply(sample_data(cm_physeq), function (y) if(class(y)!= "factor") as.factor(y) else y), stringsAsFactors=T) #convert all columns to factors
row.names(df) <- sample_names(cm_physeq)
sample_data(cm_physeq) <- sample_data(df)

#merge samples by sampling event to combine all Niskin replicates
merged <- merge_samples(cm_physeq, "Sample")
dfmerged <- merged %>% psmelt()

#check that it actually summed the taxa
ASVnames10 = names(sort(taxa_sums(cm_physeq), TRUE)[1:10]) #first 10 ASVs
CM10  = prune_taxa(ASVnames10,  cm_physeq) #phyloseq object of the first 10 ASVs and occurrence in each replicate
mCM10 = prune_taxa(ASVnames10, merged) #phyloseq object of the first 10 ASVs and summed occurrence after merging replicates
gc_samples = sample_names(subset(sample_data(cm_physeq),Sample == "GC354-527")) #take sample names from phyloseq object
print(gc_samples)
otu_table(CM10)[, gc_samples] #table of ASVs in each replicate from og phyloseq
rowSums(otu_table(CM10)[, gc_samples]) #summed ASVs from og phyloseq
otu_table(mCM10)["GC354-527", ] #table of ASV abundances from merged phyloseq
#rowSums result matches otu_table result -- merge_samples worked!

#calculate relative abundance
cm_physeq_rel <- merged %>% transform_sample_counts(function(x) {x/sum(x)})
dfrel <- cm_physeq_rel %>% psmelt()

#making sure the function worked by checking first ASV in dfrel
stetson <- dfmerged %>% subset(Sample == "Stetson-55") #subset all the ASVs from that sampling event
sum(stetson$Abundance) #total abundance from sampling event
46651/49529 #merged abundance of that first ASV / total abundance from sampling event
#matches rel abundance in dfrel -- function worked!

#WORKING ON IT -- add rel abundance to og physeq
rel_abdc <- dfrel %>% select(c(OTU, Sample, Abundance))
rel_abdc <- dplyr::rename(rel_abdc, c(sample_Sample = "Sample"))

test <- left_join(df_cm_physeq, rel_abdc, by = c("OTU", "sample_Sample"))
test$Abundance.x <- test %>% select(-c(Abundance.x))
test <- dplyr::rename(test, c(Abundance = "Abundance.y"))

exp <- otu_table(cm_physeq)

#export table of rel abundance ASVs so we can put it into Geneious to make a phylogenetic tree
tax_tab_rel <- tax_table(cm_physeq_rel) #take component of ASVs out of phyloseq object
write.csv(tax_tab_rel, "~/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_tab_rel.csv") #save as csv file

#add rel abundance tree of ctenophores and medusozoans to phyloseq object
rel_tree <- read.nexus("/Users/quattrinia/Desktop/AW RStudio/data/gomx-cteno-medu-dist/tax_tab_rel alignment FastTree Tree.nex") #import Geneious tree of aligned ctenos and medus
rooted_rel_tree <- midpoint.root(rel_tree)
rel_physeq_tree <- merge_phyloseq(cm_physeq_rel, rooted_rel_tree) #merge tree with existing phyloseq object

#calculate Unifrac distance matrix 
uni_matrix_rel <- UniFrac(rel_physeq_tree, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #unweighted UniFrac
uni_matrix_wt_rel <- UniFrac(rel_physeq_tree, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE) #weighted UniFrac

#-------NOT DOING ANYMORE plot dbRDA using gg_ordiplot, plots ordination colored by site and shaped by method-----
ordiplot <- gg_ordiplot(dbrda, groups = norm_metadata$Site, pt.size = 1.5, hull = FALSE, spiders = FALSE, ellipse = FALSE)
names(ordiplot) #list components of ordiplot object
ord.data <- ordiplot$df_ord #capture the df of point coordinates for ordination
ord.data$`CTD/ROV` <- norm_metadata$`CTD/ROV` #add sampling method to ordination df
ord.data <- dplyr::rename(ord.data, c(Site = "Group", Sampling_Method = "CTD/ROV")) #rename columns in ordination df
ggplot(data = ord.data, aes(x=x, y=y, color = Site, shape = Sampling_Method)) +
  geom_point(size = 2) +
  labs(shape = "Sampling Method", x = "dbRDA1 (15.9%)", y = "dbRDA2 (10.6%)") +
  theme_bw() 
ord.plot <- ordiplot$plot #selecting plot object 
ord.plot + theme_bw() + aes(x=x, y=y, color = Site)

#--------NOT DOING ANYMORE merging sample replicates---------

#convert sample data table columns to factor so we can merge samples
df <- as.data.frame(lapply(sample_data(cm_physeq), function (y) if(class(y)!= "factor") as.factor(y) else y), stringsAsFactors=T) #convert all columns to factors
row.names(df) <- sample_names(cm_physeq)
sample_data(cm_physeq) <- sample_data(df)

#merge samples by sampling event to combine all Niskin replicates
merged <- merge_samples(cm_physeq, "Sample")
dfmerged <- merged %>% psmelt()

#check that it actually summed the taxa
ASVnames10 = names(sort(taxa_sums(cm_physeq), TRUE)[1:10]) #first 10 ASVs
CM10  = prune_taxa(ASVnames10,  cm_physeq) #phyloseq object of the first 10 ASVs and occurrence in each replicate
mCM10 = prune_taxa(ASVnames10, merged) #phyloseq object of the first 10 ASVs and summed occurrence after merging replicates
gc_samples = sample_names(subset(sample_data(cm_physeq),Sample == "GC354-527")) #take sample names from phyloseq object
print(gc_samples)
otu_table(CM10)[, gc_samples] #table of ASVs in each replicate from og phyloseq
rowSums(otu_table(CM10)[, gc_samples]) #summed ASVs from og phyloseq
otu_table(mCM10)["GC354-527", ] #table of ASV abundances from merged phyloseq
#rowSums result matches otu_table result -- merge_samples worked!

#get merged otu table back into phyloseq format
merge_otu <- as.data.frame(otu_table(merged))
transpose <- data.frame(t(merge_otu)) #switch row and column names
transpose <- dplyr::rename(transpose, c("Bright-111" = "Bright.111", "Bright-67" = "Bright.67", "Bright-84-Background" = "Bright.84.Background",
                                        "Bright-85-Swiftia" = "Bright.85.Swiftia", "Bright-CTD-1.4" = "Bright.CTD.1.4", "Bright-CTD-23.8" = "Bright.CTD.23.8",
                                        "Bright-CTD-44.7" = "Bright.CTD.44.7", "Bright-CTD-5.8" = "Bright.CTD.5.8", "Bright-CTD-77.45" = "Bright.CTD.77.45",
                                        "GC354-527" = "GC354.527", "GC354-531" = "GC354.531", "Stetson-55" = "Stetson.55", "VK826-474" = "VK826.474",
                                        "VK826-CTD-1.058" = "VK826.CTD.1.058", "VK826-CTD-20.961" = "VK826.CTD.20.961", "VK826-CTD-5.214" = "VK826.CTD.5.214",
                                        "VK826-CTD-9.577" = "VK826.CTD.9.577")) #rename column names because got changed from dashes to decimals
otu_table(merged) <- otu_table(transpose, taxa_are_rows = TRUE) #put back into phyloseq object

#put cm_physeq sample data table into merged phyloseq, because merging makes the sample data weird 
og_sample <- as.data.frame(sample_data(cm_physeq))
#since most sample data is the same across replicates (except pH), taking the first replicate as a representative of the sampling data for all replicates that were merged
rep1 <- og_sample[c("trimmed.1368.01.SH18291.meso.octo_S1_L001_R1_001.fastq.gz", "trimmed.1368.05.SH18295.meso.octo_S5_L001_R1_001.fastq.gz", "trimmed.1368.08.SH18299.meso.octo_S8_L001_R1_001.fastq.gz","trimmed.1368.11.SH18302.meso.octo_S11_L001_R1_001.fastq.gz",
                    "trimmed.1368.14.SH18305.meso.octo_S14_L001_R1_001.fastq.gz","trimmed.1368.16.SH18307.meso.octo_S16_L001_R1_001.fastq.gz","trimmed.1368.18.SH18309.meso.octo_S18_L001_R1_001.fastq.gz","trimmed.1368.21.SH18312.meso.octo_S21_L001_R1_001.fastq.gz",
                    "trimmed.1368.26.SH18353.meso.octo_S26_L001_R1_001.fastq.gz","trimmed.1368.29.SH18357.meso.octo_S29_L001_R1_001.fastq.gz","trimmed.1368.33.SH18180.deep.octo_S33_L001_R1_001.fastq.gz","trimmed.1368.36.SH18196.deep.octo_S36_L001_R1_001.fastq.gz",
                    "trimmed.1368.39.SH18199.deep.octo_S39_L001_R1_001.fastq.gz","trimmed.1368.42.SH18202.deep.octo_S42_L001_R1_001.fastq.gz","trimmed.1368.44.SH18204.deep.octo_S44_L001_R1_001.fastq.gz","trimmed.1368.47.SH18376.deep.octo_S47_L001_R1_001.fastq.gz",
                    "trimmed.1368.51.SH18380.deep.octo_S51_L001_R1_001.fastq.gz"), ] #subset first replicate from each sampling event
row.names(rep1) <- rep1$Sample #make samples the row names to match otu table column names
sample_data(merged) <- sample_data(rep1) #put back into merged phyloseq

#MIGHT CHANGE -- select only one column from each sampling event to match dimensions of distance matrix, since variable values are mostly the same
norm_metadata$Row <- row.names(norm_metadata)
norm_metadata <- norm_metadata %>% subset(Row == "1"|Row == "5"|Row == "8"|Row == "11"|Row == "14"|Row == "16"|Row == "18"|Row == "20"|Row == "24"
                                          |Row == "27"|Row == "29"|Row == "32"|Row == "35"|Row == "38"|Row == "40"|Row == "42"|Row == "46")


