# Load and install packages
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("phyloseq", "ComplexHeatmap"), update = FALSE)

library("phyloseq")

setwd("C:/Users/Marcelo/Desktop/EMP_soil_strain_analysis")

##### Load metadata files

# Get list of metadata files
metadata_files <- list.files("./metadata")

metadata_files

# Join all metadata files into a single DF
# iterate over all files
for(file_number in seq(from = 1, to = length(metadata_files), by = 1)){
  # If it is the first file, use it create the result DF.
  if(file_number == 1){
    sample_metadata <- readr::read_tsv(paste("./metadata", metadata_files[file_number],
                                      sep = "/"))
    print(ncol(sample_metadata))
    print(nrow(sample_metadata))
  }
  # rbind the rest of files
  else{
    current_table <- readr::read_tsv(paste("./metadata", metadata_files[file_number],
                                    sep = "/"))
    sample_metadata <- rbind(sample_metadata, current_table)
  }
}

head(sample_metadata)

# Number of unique variables (columns)
length(unique(sample_metadata))

# Number of unique samples (rows)
length(unique(sample_metadata$sample_name))

write.csv(sample_metadata, "sample_metadata.csv", row.names = FALSE, quote = FALSE)

##### Let´s filter samples

# Filter only Soil samples (empo_4 == "Soil (non-saline)") and biomes with anthropogenic influence (env_biome).
anthropogenic_biomes <- c("rangeland biome", "anthropogenic terrestrial biome", "urban biome", "cropland biome")

soil_metadata <- sample_metadata %>%
  filter(empo_4 == "Soil (non-saline)") %>%
  filter(env_biome %in% anthropogenic_biomes)

head(soil_metadata)

# Number of unique variables (columns)
length(unique(soil_metadata))

# Number of unique samples (rows)
length(unique(soil_metadata$sample_name))

write.csv(sample_metadata, "soil_anthro_metadata.csv", row.names = FALSE, quote = FALSE)

# Get the names only of the samples that we choose
soil_sample_names <- soil_metadata$sample_name

##### Read OTU tables

# get list of biom files
biom_files <- list.files("./biom_files")

biom_files

for (i in seq(from = 1, to = length(biom_files), by = 1)) {
  print(biom_files[i])
  x <- biom_files[i]
  biom_path <- paste("./biom_files", x, sep = "/")
  print(biom_path)
  eval(call("<-", as.name(x), phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_greengenes)))
}

# Merge all biom files into a single one.
biom_merged <- do.call(merge_phyloseq, mget(biom_files))

biom_merged

# Collapse samples at Genus level
biom_genus <- phyloseq::tax_glom(biom_merged, "Genus", NArm = TRUE)

# Extract otu table and tax table from merged biom file and cbind them
asv_table <- cbind(data.frame(otu_table(biom_genus)), data.frame(tax_table(biom_genus)))

# Lets change the name of the samples in original metadata file because they don't comply with R standards.
soil_sample_names <- paste0("X", soil_sample_names)

# Lets filter the otu_table with taxonomy and choose only the selected samples.
soil_table <- asv_table %>% select(any_of(c(soil_sample_names, c("Kingdom", "Phylum",  "Class",  "Order",  "Family", "Genus", "Species"))))

soil_table <- soil_table %>% 
  unite(bacteria, c("Kingdom", "Phylum",  "Class",  "Order",  "Family", "Genus", "Species"), sep = "_")

# moving tax column to the first column
soil_table <- cbind(soil_table[, ncol(soil_table)], soil_table[1:nrow(soil_table), 1:(ncol(soil_table) - 1)])

# renaming tax to taxonomy. rename() is a dplyr function.
soil_table <- dplyr::rename(soil_table, taxonomy = "soil_table[, ncol(soil_table)]")

# Collapse all ASV entries with the same taxonomic id
soil_table <- soil_table %>% group_by(taxonomy) %>% summarise_all(list(sum))

# setting row names and dropping rownames column
soil_table <- tibble::remove_rownames(soil_table)
soil_table <- tibble::column_to_rownames(soil_table, var = "taxonomy")

head(soil_table)

soil_table_ordered_by_means <- soil_table[order(rowMeans(soil_table), decreasing = TRUE),]

# Select only the n more abundant species.
soil_table_ordered_50 <- soil_table_ordered_by_means[1:50,]

# Remove columns (samples) with yero count
soil_table_ordered_50 <- soil_table_ordered_50[, colSums(soil_table_ordered_50 != 0) > 0]

# Generate a column with the names of ASVs/OTUs using rownames.
soil_table_bp <- soil_table_ordered_50
soil_table_bp["bacteria"] <- row.names(soil_table_bp)

colnames(soil_table_bp)

# Gather
soil_table_bp <- gather(soil_table_bp, X13114.control.soil.grass.near.BRF:X13114.shade.23.s010 , key = "sample", value = "abundance")

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

ggplot(soil_table_bp, aes(x=sample, y=abundance, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=cbbPalette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

############ Heatmaps #

library(ComplexHeatmap)

#scaled by column (samples)
table_hm_c_scaled <- scale(soil_table_ordered_50)
#scaled by row (species)
table_hm_r_scaled <- scale(t(soil_table_ordered_50))

# Graph heatmap samples-scaled
heatmap(table_hm_c_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")
# Graph heatmap species-scaled
heatmap(table_hm_r_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

heatmap(data.matrix(soil_table_ordered_50), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale = "column")

### Heatmap with annotations

biomes_metadata <- select(soil_metadata, sample_name, env_biome)

biomes_metadata$sample_name = paste0("X", biomes_metadata$sample_name)

biomes_metadata <- biomes_metadata %>%
  filter(sample_name %in% colnames(soil_table_ordered_50))

biomes <- biomes_metadata$env_biome

column_ha = ComplexHeatmap::HeatmapAnnotation(biome = biomes)

# 
ComplexHeatmap::Heatmap(table_hm_c_scaled, name = "Soil genus-level ASVs abundance", top_annotation = column_ha)

# 
ComplexHeatmap::Heatmap(t(table_hm_r_scaled), name = "Soil genus-level ASVs abundance", top_annotation = column_ha)
