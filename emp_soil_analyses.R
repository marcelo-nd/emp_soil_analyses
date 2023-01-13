# Load and install packages
library(dplyr, quietly = TRUE)
library(readr)
library(tidyr)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("phyloseq", quietly = TRUE))
  BiocManager::install("phyloseq", update = FALSE)

if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap", update = FALSE)

library("phyloseq", quietly = TRUE)

setwd("C:/Users/Marcelo/OneDrive - UT Cloud/Postdoc Tü/Sci/EMP_soil_strain_analysis")

source("C:/Users/Marcelo/Documents/Github/microbiome-help/table_importers.R")
source("C:/Users/Marcelo/Documents/Github/microbiome-help/graphs.R")

##################################################################################
# Processing data files from the EMP

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
                                      sep = "/"), name_repair = "unique", col_types = cols())
  }
  # rbind the rest of files
  else{
    current_table <- readr::read_tsv(paste("./metadata", metadata_files[file_number],
                                           sep = "/"), name_repair = "unique", col_types = cols())
    if (identical(colnames(sample_metadata), colnames(current_table))) {
      sample_metadata <- rbind(sample_metadata, current_table)
    }
  }
}

sample_metadata$sample_name <- make.names(sample_metadata$sample_name)

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

write.csv(soil_metadata, "soil_anthro_metadata.csv", row.names = FALSE, quote = FALSE)

# Get the names only of the samples that we choose
soil_sample_names <- soil_metadata$sample_name

#soil_sample_names <- make.names(soil_sample_names)

soil_sample_names

##### Read OTU tables

# get list of biom files
biom_files <- list.files("./biom_files/150")

biom_files

# Import all biom files in directory.
for (i in seq(from = 1, to = length(biom_files), by = 1)) {
  print(biom_files[i])
  biom_file_name <- biom_files[i]
  biom_path <- paste("./biom_files/150", biom_file_name, sep = "/")
  print(biom_path)
  eval(call("<-", as.name(biom_file_name), phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_greengenes)))
}

# Merge all biom files into a single one.
biom_merged <- do.call(merge_phyloseq, mget(biom_files))

biom_merged


##################################################################################
##### Genus Level

# Extract and clean ASV table.

asv_table <- extract_table(biom_object = biom_merged, tax_rank = "Genus", col_names = get_colNames_per_rank("Genus"))

asv_table <- clean_table(feature_table = asv_table, order_table = TRUE)

# Select data only from soil samples.
soil_table <- asv_table %>% select(any_of(soil_sample_names))

head(soil_table)

# Select only the n more abundant ASVs.
soil_table_ordered_50 <- soil_table[1:50,]

# Remove columns (samples) with zero count
soil_table_ordered_50 <- soil_table_ordered_50[, colSums(soil_table_ordered_50 != 0) > 0]

# Export the table
write.table(soil_table_ordered_50, "soil_asv_table.csv", row.names = TRUE, quote = FALSE, col.names = FALSE, sep = ",")

# Do barplot.
barplot_from_feature_table(soil_table_ordered_50)

##################################################################################
#Heatmaps

library(ComplexHeatmap)

# Scaled by column (samples)
table_hm_c_scaled <- scale(soil_table_ordered_50)
# Scaled by row (species)
table_hm_r_scaled <- t(scale(t(soil_table_ordered_50)))

### Heatmap with annotations
# Select sample_name and env_biome columns from metadata for each sample.
biomes_metadata <- select(soil_metadata, sample_name, env_biome)

# Correct samples names
#biomes_metadata$sample_name = paste0("X", biomes_metadata$sample_name)

# Get samples from metadata that are also in samples used to barplot (samples with "counts = 0" were removed)
biomes_metadata <- biomes_metadata %>%
  filter(sample_name %in% colnames(soil_table_ordered_50))

head(biomes_metadata)

# Check that names in scaled matrix and metadata are identical so that metadata info matches the samples in abundance matrix.
identical(colnames(table_hm_c_scaled), biomes_metadata$sample_name)

# Create annotation for heatmap using env_biome orign.
column_ha = ComplexHeatmap::HeatmapAnnotation(biome = biomes_metadata$env_biome)

# Create function name
col_fun = circlize::colorRamp2(c(min(table_hm_c_scaled), 0.5, max(table_hm_c_scaled)), c("#5F8D4E", "white", "#FF6464"))

# Heatmap, column-scaled
ComplexHeatmap::Heatmap(table_hm_c_scaled,
                        name = "Soil genus-level ASVs abundance",
                        top_annotation = column_ha,
                        show_column_names = FALSE,
                        col = col_fun)

# Heatmap, row-scaled Remove rows (ASVs) with sum counts of 0
table_hm_r_scaled <- table_hm_r_scaled[apply(table_hm_r_scaled[,-1], 1, function(x) !all(is.nan(x))),]
ComplexHeatmap::Heatmap(table_hm_r_scaled,
                        name = "Soil genus-level ASVs abundance",
                        top_annotation = column_ha,
                        show_column_names = FALSE)

##################################################################################
# Co-ocurrence network
library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(soil_table_ordered_50, method = "pa")

# Infering co-ocurrences
cooccur.otus <- cooccur::cooccur(otu_table_pa,
                        type = "spp_site",
                        spp_names = TRUE)

summary(cooccur.otus)
plot(cooccur.otus)


# Getting only the significant interactions
co <- print(cooccur.otus)

# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(otu_table_pa),
                    label = rownames(otu_table_pa),
                    color = "#606482",
                    shadow = TRUE) 

# Create an edges dataframe from the significant pairwise co-occurrences.
edges <- data.frame(from = co$sp1, to = co$sp2,
                    color = ifelse(co$p_lt <= 0.05,
                                   "red", "green"),
                                   dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

# Plotting network
library(visNetwork)
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")

g <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
print(g, e=TRUE, v=TRUE)
plot(g)
write.graph(g, "C:/Users/marce/Desktop/coocur_network.graphml", format = "graphml")
##################################################################################
##### Family Level
# Extract ASV table

asv_table_fam <- extract_table(biom_object = biom_merged, tax_rank = "Family", col_names = get_colNames_per_rank("Family"))

asv_table_fam <- clean_table(feature_table = asv_table_fam, order_table = TRUE)

soil_table_fam <- asv_table_fam %>% select(any_of(soil_sample_names))

head(soil_table_fam)

# Select only the n more abundant ASVs.
soil_table_50_fam <- soil_table_fam[1:50,]

barplot_from_feature_table(soil_table_50_fam)

##################################################################################

