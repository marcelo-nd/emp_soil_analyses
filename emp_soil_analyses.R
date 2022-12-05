# Load and install packages

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

library("phyloseq")

##### Load metadata files

# Get list of metadata files
metadata_files <- list.files("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis/metadata")

# Iterate over each file.
for(file_number in seq(from = 1, to = length(metadata_files), by = 1)){
  
  if(file_number == 1){
    sample_metadata <- read_delim(paste("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis/metadata", metadata_files[file_number], sep = "/"),
                                                    delim = "\t", escape_double = FALSE,
                                                    trim_ws = TRUE)
  }
  else{
    sample_metadata <- rbind(sample_metadata, read_delim(paste("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis/metadata", metadata_files[file_number], sep = "/"),
                                                        delim = "\t", escape_double = FALSE,
                                                        trim_ws = TRUE))
  }
}

head(sample_metadata)

soil_metadata <- sample_metadata %>% filter(empo_4 == "Soil (non-saline)")

# Selecting only samples from biomes with anthropogenic influence.
biomes <- c("rangeland biome", "anthropogenic terrestrial biome", "urban biome", "cropland biome")
soil_metadata <- soil_metadata %>% filter(env_biome %in% biomes)


head(soil_metadata)

soil_sample_names <- soil_metadata$sample_name


# Now read otu tables

biom_files <- list.files("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis/biom_files")

biom_files

files_list <- c()

for (i in seq(from = 1, to = length(biom_files), by = 1)) {
  print(biom_files[i])
  x <- biom_files[i]
  biom_path <- paste("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis/biom_files", x, sep = "/")
  print(biom_path)
  eval(call("<-", as.name(x), import_biom(biom_path, parseFunction=parse_taxonomy_greengenes)))
  files_list <- c(files_list, biom_files[i])
}

files_list

biom_merged <- do.call(merge_phyloseq, mget(files_list))

otu_table <- cbind(data.frame(otu_table(biom_merged)), data.frame(tax_table(biom_merged)))

soil_sample_names2 <- c()

for (i in seq(from = 1, to = length(soil_sample_names), by = 1)){
  soil_sample_names2 <- c(soil_sample_names2, paste("X", soil_sample_names[i], sep = ""))
}

soil_sample_names2

otu_table <- otu_table %>% select(any_of(c(soil_sample_names2, c("Kingdom", "Phylum",  "Class",  "Order",  "Family", "Genus", "Species"))))

head(otu_table)

otu_table <- otu_table %>% 
  unite(bacteria, c("Kingdom", "Phylum",  "Class",  "Order",  "Family", "Genus", "Species"))

head(otu_table)

# moving tax column to the first column
otu_table2 <- cbind(otu_table[, ncol(otu_table)], otu_table[1:nrow(otu_table), 1:(ncol(otu_table) - 1)])

head(otu_table2)

# renaming tax to taxonomy. rename() is a dplyr function.
otu_table2 <- dplyr::rename(otu_table2, taxonomy = "otu_table[, ncol(otu_table)]")

otu_table2["taxonomy"] <- dereplicate_taxonomy(otu_table2$taxonomy)

# setting row names and dropping rownames column
otu_table2 <- tibble::remove_rownames(otu_table2)
otu_table2 <- tibble::column_to_rownames(otu_table2, var = "taxonomy")

head(otu_table2)

otu_table_ordered_means <- otu_table2[order(rowMeans(otu_table2), decreasing = TRUE),]

# Select only the 30 more abundant species.
otu_table3 <- otu_table_ordered_means[1:30,]

# Generate a column with the names of ASVs/OTUs using rownames.
otu_table3["bacteria"] <- row.names(otu_table3)

colnames(otu_table3,)

# For all soil samples
otu_g <- gather(otu_table3, X13114.king.27.s001:X13114.stegen.38.s018 , key = "sample", value = "counts")

# For soil samples with anthopogenic influence
otu_g <- gather(otu_table3, X13114.control.soil.grass.near.BRF:X13114.shade.23.s010 , key = "sample", value = "counts")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2","brown1", "#CC79A7", "olivedrab3", "rosybrown",
                "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
                "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
                "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
                "#e3ae78", "#a23f3f", "#290f76", "#ce7e00", "#386857",
                "#738564", "#e89d56", "#cd541d", "#1a3a46", "#ffe599",
                "#583E26", "#A78B71", "#F7C815", "#EC9704", "#9C4A1A",
                "firebrick2", "#C8D2D1", "#14471E", "#EE9B01", "#DA6A00",
                "#4B1E19", "#C0587E", "#FC8B5E", "#EA592A", "#FEF4C0")

ggplot(otu_g, aes(x=sample, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)

+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#scaled by column
otu_table_scaled <- scale(otu_table_ordered_means[1:30,])

# Graph heatmap
heatmap(otu_table_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

species_scaled_df <- scale(t(otu_table_ordered_means[1:30,]))
heatmap(t(species_scaled_df), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

heatmap(data.matrix(otu_table2), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale = "column")


unique(soil_metadata$env_biome)
