# HMP BacDive Trait Extraction Script
# Extracts bacterial traits from BacDive database for Human Microbiome Project species
# Author: Feargal Ryan (https://github.com/feargalr)
# Date: August 19, 2025

library(BacDive)
library(curatedMetagenomicData)
source("Helper Functions.R")

# Extract species names to search BacDive using Human Microbiome Project dataset
# Load HMP 2012 relative abundance data from curatedMetagenomicData
HMP_data = curatedMetagenomicData(pattern = "2021-03-31.HMP_2012.relative_abundance",counts = TRUE,dryrun = FALSE)
HMP_data = HMP_data$`2021-03-31.HMP_2012.relative_abundance`

# Extract rownames which contain taxonomic information
species_names = rownames(HMP_data)
# Filter for species-level annotations (s__ prefix)
species_names = species_names[grepl("s__",species_names)]
species_names = unique(species_names)

# Parse taxonomic strings to extract species names
# Taxonomic format: Kingdom|Phylum|Class|Order|Family|Genus|Species
species_names = sapply(as.character(species_names),function(y) {strsplit(x = y,split="\\|")[[1]][7]})
# Remove species prefix "s__"
species_names = gsub("s__","",species_names)
names(species_names) = NULL
species_names <- unique(species_names)

# Clean up species names for BacDive queries
# Remove metagenome-assembled genome identifiers and species placeholders
species_names = species_names[!grepl("_CAG_",species_names)]  # Remove CAG (Co-Abundance Group) entries
species_names = species_names[!grepl("_sp_",species_names)]   # Remove generic "sp." entries
# Convert underscores to spaces for proper species names
species_names = gsub("_"," ",species_names)

# Initialize counter for successful Gram stain retrievals
counter = 0

# Initialize BacDive API connection
# Sign up for API access at https://api.bacdive.dsmz.de
bd <- open_bacdive(
  username = "",  # Add your BacDive username here
  password = ""   # Add your BacDive password here
)

# Initialize collector for all trait entries across species in long format
trait_long_df <- data.frame()

# Main processing loop - iterate through each species
for (sp in species_names) {
  message("Processing: ", sp)
  
  # Attempt to retrieve records from BacDive with error handling
  records <- tryCatch(retrieve(bd, query = sp, sleep = 0.1), error = function(e) NULL)
  
  # Skip species if no records found
  if (is.null(records) || length(records) == 0) {
    cat("Species:", sp, "\n  ❌ No records found\n\n")
    next
  }
  
  # Select most authoritative record (preferring type strains)
  idx <- get_type_strain_index(records)
  preferred <- get_preferred_record(records)
  
  # Extract strain designation information
  strain <- preferred$`Name and taxonomic classification`$`strain designation`
  strain_str <- if (is.null(strain) || all(trimws(unlist(strain)) == "")) {
    "not reported"
  } else {
    paste(unlist(strain), collapse = "; ")
  }
  
  # Extract morphological traits using helper functions
  gram_stain <- extract_morphology_trait(preferred, "gram stain")
  cell_shape <- extract_morphology_trait(preferred, "cell shape")
  motility   <- extract_morphology_trait(preferred, "motility")
  
  # Count successful Gram stain extractions for success rate calculation
  if(gram_stain != "not reported"){counter = counter+1}
  
  # Extract physiological traits
  metabolites <- extract_positive_metabolites(preferred)  # Metabolites the organism can utilize
  enzymes     <- extract_positive_enzymes(preferred)      # Enzymes with positive activity
  oxygen_tolerance <- extract_oxygen_tolerance(preferred) # Aerobic/anaerobic classification
  
  # Optional output for debugging (currently commented out)
  #cat("  Positive metabolites: ", paste(metabolites, collapse = ", "), "\n")
  #cat("  Positive enzymes:     ", paste(enzymes, collapse = ", "), "\n\n")
  #cat("  Oxygen tolerance:  ", oxygen_tolerance, "\n")
  #cat("Species:", sp, "\n")
  #cat("  Type strain found? ", ifelse(length(idx) > 0, "✅ YES", "❌ NO"), "\n")
  #cat("  Selected strain:   ", strain_str, "\n")
  #cat("  Gram stain:        ", gram_stain, "\n")
  #cat("  Cell shape:        ", cell_shape, "\n")
  #cat("  Motility:          ", motility, "\n\n")
  
  # Build long-format data entries for this species
  # Each trait type gets separate rows for downstream analysis
  rows <- list()
  
  # Add trait data in long format (one row per trait value)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "metabolite", metabolites)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "enzyme", enzymes)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "oxygen", "oxygen_tolerance", oxygen_tolerance)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "gram_stain", "gram_stain", gram_stain)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "cell_shape", "cell_shape", cell_shape)
  rows[[length(rows)+1]] <- append_trait_rows(sp, "motility", "motility", motility)
  
  # Append to main data frame
  trait_long_df <- rbind(trait_long_df, do.call(rbind, rows))
}

# Calculate and display success rate for Gram stain extraction
result = round((counter/length(species_names))*100,2)
result = paste(result,"%",sep="")
textfm = "Found result for "
print(paste(textfm,result))
