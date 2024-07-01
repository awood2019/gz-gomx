#comprehensive classifier of Cnidaria and Ctenophora in Gulf of Mexico 28S eDNA

#--------code from GoMx phyla distribution script to make phyloseq object combining taxonomy, sample metadata, and abundance data-------

#load packages
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 
library("ggplot2")

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-phy-dist")

#load components of phyloseq object: taxonomy table, count table, and sample data table.
tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & phyla
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_phy: includes taxonomic information for each representative (ASV) sequence
phylum <- tax_tab$Phylum #pulling out phylum column from taxonomy table
tax_tab_phy <- tibble(phylum) #making phyla into a tibble containing phylum for each sequence
tax_tab_phy <- as.matrix(tax_tab_phy) #make tibble into matrix
row.names(tax_tab_phy) <- tax_tab$Sequence #make sequence column the row names
tax_tab_phy[is.na(tax_tab_phy)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_phy: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_phy)
count_tab_phy <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_phy) <- count_tab$...1 #make sequences the row names (ignore warning message)

#sample_info_tab_phy: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab_phy <- sample_info_tab
sample_info_tab_phy <- sample_info_tab_phy[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_phy) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object with just the count table and taxonomy table
ASV_physeq <- phyloseq(otu_table(count_tab_phy, taxa_are_rows = TRUE), tax_table(tax_tab_phy), sample_data)
ASV_physeq <- prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'allphy_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe for easy viewing
df_ASV_physeq <- ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe
head(df_ASV_physeq)

#----------end of the code from GoMx phyla distribution script--------

#---------assign taxonomy to Cnidaria and Ctenophora ASVs--------

#load packages
library(dada2)

#set working directory
setwd("~/Desktop/AW RStudio/data/gomx-cnid-cteno-classifier")

#select only the cnidarians and ctenophores
cnid_cteno_table <- df_ASV_physeq %>%
  subset(phylum == "Cnidaria" | phylum == "Ctenophora") %>%
  filter(Abundance > 0)
write.table(cnid_cteno_table, file = "Cnid_cteno_preclassifier.tsv", sep = "\t",row.names = FALSE, quote = FALSE) #save file as a tsv

#create vector of only unique sequences from table so we can compare them to the FASTA files
seqs <- unique(cnid_cteno_table$OTU)

#assignTaxonomy using FASTA file as reference
taxa <- assignTaxonomy(seqs, "28S_Cnid_Cteno_assignTaxonomy.fasta", multi=TRUE, minBoot = 80,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) #classifying to the genus level with assignTaxonomy
unname(taxa)
taxa <- taxa[,-c(7)] #remove empty column
unique(taxa [,6]) #see which ASVs were identified to genus level

#assignSpecies using FASTA file as reference
species <- addSpecies(taxa, "28S_Cnid_Cteno_assignSpecies.fasta", allowMultiple=TRUE) #finding 100% matches to our reference database of Gulf of Mexico ctenos with assignSpecies
unique(species [, 7]) #see how many were identified to species level
species[species == ''] <- NA

#convert from vector to dataframe
species_df <- as.data.frame(species)
species_df$seq <- row.names(species_df)

#add abundance of eaech unique ASV to dataframe
length(unique(species_df[["seq"]])) #find how many unique ASVs are in table

#make table of abundance counts for each unique ASV
counts_per_ASV <- cnid_cteno_table %>%
  group_by(OTU) %>%
  summarize(totalcount = sum(Abundance))
names(counts_per_ASV)[names(counts_per_ASV) == "OTU"] <- "seq"

#combine abundance table with classifier table
species_df <- left_join(species_df, counts_per_ASV, by = "seq")

#save as a table
setwd("~/Desktop/AW RStudio/results/gomx-cnid-cteno-classifier")
write.table(species_df, file = "Cnid_Cteno_classifier_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#create table from metadata table with only depth, method, and abundance
cnid_cteno_depth <- cnid_cteno_table %>%
  select(OTU, Abundance, depth_bin, CTD.ROV, Depth, Site)
cnid_cteno_spec_depth <- right_join(species_df, cnid_cteno_depth, by =c("seq" = "OTU")) #join with table of cnidarian and ctenophore species taxonomy
cnid_cteno_spec_depth <- cnid_cteno_spec_depth %>%    #remove negative controls
  filter(CTD.ROV !="Negative" & CTD.ROV != "NTC")
unique(cnid_cteno_depth$depth_bin)

#plot by class
ggplot(cnid_cteno_spec_depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Class)) +     #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by order
ggplot(cnid_cteno_spec_depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Order)) +     #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 



