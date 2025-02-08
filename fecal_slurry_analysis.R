library(phyloseq)
library(ggplot2)
library(microbiome)
library(mia)


### separately transformed ####

library(compositions)

otu <- read.table("feature-table.txt", sep = '\t', header = TRUE, row.names = 1,
                  check.names = FALSE)

otu_transformed <- otu + 0.1
clr_transformed <- clr(otu_transformed)
View(clr_transformed)

write.table(as.data.frame(as.matrix(clr_transformed)), file = "clr_otu_table.tsv", sep = "\t", quote = FALSE, col.names = NA)


### mia transformation from biom ####

#read and convert biom to summarized experiment (SE)
biom_object <- biomformat::read_biom("feature_table_hdf5.biom")
tse <- makeTreeSEFromBiom(biom_object)

#relative abundance
otu_rel_abd <- transformAssay(tse, method = "relabundance")

#clr transform
otu_clr <- transformAssay(tse, assay.type = "counts", method = "clr", pseudocount = TRUE)

#rclr transformation
otu_rclr <- transformAssay(tse, assay.type = "counts", method = "rclr", pseudocount = TRUE)

#storing output in variable
otu_rel_abd <- assay(otu_rel_abd, "relabundance")
otu_clr <- assay(otu_clr, "clr")
otu_rclr <- assay(otu_rclr, "rclr")

### phyloseq object ####

#reading all the tables, column 1 becomes the row names
# otu <- read.table("feature-table.txt", sep = '\t', header = TRUE, row.names = 1, 
#                    check.names = FALSE)
meta <- read.table("metadata.txt", sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
tax <- read.table("taxonomy_gg.txt", sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

#creating phyloseq object
physeq <- phyloseq((otu_table(as.matrix(otu_clr), taxa_are_rows = TRUE)), tax_table(as.matrix(tax)), sample_data(meta))

# #phyloseq object on separately clr transformed OTU table (not mia)
# physeq_clr <- phyloseq((otu_table(as.matrix(otu_clr), taxa_are_rows = TRUE)), tax_table(as.matrix(tax)), sample_data(meta))


### abundance-based filtering ####

#abundance threshold (1% relative abundance)
abundance_threshold <- 0.1

#relative abundance of each taxon across all samples
taxa_abundance <- apply(otu_table(physeq), 1, mean)

#filter taxa based on the abundance threshold
abundant_taxa <- names(taxa_abundance[taxa_abundance >= abundance_threshold])

#phyloseq object to include only the abundant taxa
physeq_abundant <- prune_taxa(abundant_taxa, physeq)

cat("Number of taxa retained after abundance-based filtering:", ntaxa(physeq_abundant), "\n")


### stacked bar plots - phylum ####

#removing empty phylum & unassigned domain
phy1 <- subset_taxa(physeq, !(Domain %in% c("Unassigned")))
phy1 <- subset_taxa(phy1, !(Phylum %in% c("p__")))
 
#subset sample based on a certain criteria of a column
phy1 <- subset_samples(phy1, replicate == "replicate_2")

#calculating relative abundance
phy1_rel_abd <- transform_sample_counts(phy1, function(x) (x/sum(x)))

#normal plot                    
plot_bar(phy1, fill = "Phylum") + 
  geom_bar(aes(color = Phylum), 
           stat = "identity", position = "stack")

#phylum in different groups and x-axis changed to the names in column "Sample_name"
plot_bar(phy1, fill = "Phylum") + 
  aes(x = Sample_name, color = Phylum) + 
  facet_grid(~donor, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "stack")


### stacked bar plots - genus ####

#subset sample based on a certain criteria of a column & taxa
phy2 <- subset_taxa(physeq, !(Genus %in% c("g__")))
phy2 <- subset_samples(phy2, replicate == "replicate_1")

phy2 <- tax_glom(phy2, "Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
phy2 <- subset_taxa(phy2, !is.na(Genus) & Genus != "")

#calculate the total abundance of each genus across all samples
genus_abundance <- apply(otu_table(phy2), 1, median)

#get the top 20 genera based on total abundance
top_20_genera <- names(sort(genus_abundance, decreasing = TRUE))[1:20]
physeq_top_20 <- prune_taxa(top_20_genera, phy2)

#bar plot
plot_bar(physeq_top_20, fill = "Genus") + 
  aes(x = Sample_name, color = Genus) + 
  facet_grid(~pairs, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "stack")


### normalizing w.r.t blank for each donor ####

#subset sample based on a certain criteria of a column & taxa
phy3 <- subset_taxa(physeq_abundant, !(Genus %in% c("g__")))
phy3 <- subset_samples(phy3, replicate == "replicate_1")

phy3 <- tax_glom(phy3, "Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
phy3 <- subset_taxa(phy3, !is.na(Genus) & Genus != "")

#calculate the total abundance of each genus across all samples
genus_abundance <- apply(otu_table(phy3), 1, sum)

#get the top 20 genera based on total abundance
top_20_genera <- names(sort(genus_abundance, decreasing = TRUE))[1:20]
physeq_top_20 <- prune_taxa(top_20_genera, phy3)

#separating donors
phy3a <- subset_samples(physeq_top_20, donor == "1")
phy3b <- subset_samples(physeq_top_20, donor == "2")

#sample IDs for blank
base_sample_id1 <- "FS28"
base_sample_id2 <- "FS16"
# base_sample_id1 <- "FS34"
# base_sample_id2 <- "FS31"

#convert phyloseq object into a data frame
otu_df1 <- as.data.frame(otu_table(phy3a))
otu_df2 <- as.data.frame(otu_table(phy3b))

#normalizing (ratio) w.r.t. blank
#only 1 for loop used since rows are same for both data frames
for (i in 1:nrow(otu_df1)) {
  base_value1 <- otu_df1[i, base_sample_id1]
  base_value2 <- otu_df2[i, base_sample_id2]
  otu_df1[i, -1] <- otu_df1[i, -1] / base_value1
  otu_df2[i, -1] <- otu_df2[i, -1] / base_value2
}

print(otu_df1)
print(otu_df2)

#combining data frames based on similar row names (without duplicates)
#this adds default numeric sequence (row names)
otu_df <- merge(otu_df1, otu_df2, by = "row.names")

#removing the default numeric sequence
#changing the "Row.names" column to row names then deleting the column
row.names(otu_df) <- otu_df$Row.names
otu_df$Row.names <- NULL

#creating phyloseq object with normalized data
phy3 <- phyloseq((otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)), 
                 tax_table(as.matrix(tax)), sample_data(meta))


### line plots for top genus ####

# # Select the genus of interest
# selected_genus <- "g__Prevotella"

#folder for all plots
dir.create("line_plots")
setwd("line_plots/")

#phyloseq object to data frame
df2 <- psmelt(phy3)

#creating all 20 plots
for (taxon in unique(df2$Genus)) {
  
  #declaring filename
  filename <- paste(taxon, ".png", sep = "")
  png(filename, width = 800, height = 600)
  
  #testing
  print(taxon)
  print(filename)
  
  #subset the phyloseq object to include only the selected genus
  physeq_selected <- subset_taxa(phy3, Genus == taxon)
  
  #convert the phyloseq object to a data frame
  df <- psmelt(physeq_selected)
  
  #line and point graph
  plot <- ggplot(df, aes(x = Sample_name, y = Abundance, group = 1)) +
    geom_line() +
    geom_point() +
    labs(title = paste("Abundance of", taxon, "across samples and donors"),
         x = "Sample", y = "Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~donor, scales = "free", space = "free")
  
  #save the image
  print(plot)
  dev.off()
}

setwd("~/debaleena/collabs/SRC_fecal_slurry/")

