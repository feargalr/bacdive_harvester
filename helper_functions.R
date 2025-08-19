# BacDive Helper Functions
# Author: Feargal Ryan (https://github.com/feargalr)
# Date: August 19, 2025

# Find indices of type strain records in a list of bacterial records
# Type strains are reference strains that define a bacterial species
get_type_strain_index <- function(record_list) {
  which(sapply(record_list, function(rec) {
    # Navigate to the type strain field in the taxonomic classification
    val <- rec$`Name and taxonomic classification`$`type strain`
    # Check if the value exists and equals "yes" (case-insensitive)
    !is.null(val) && tolower(val) == "yes"
  }))
}

# Select the most authoritative record from a list of bacterial records
# Prioritizes type strains, falls back to first record if no type strain exists
get_preferred_record <- function(record_list) {
  idx <- get_type_strain_index(record_list)
  if (length(idx) > 0) {
    # Return first type strain found
    return(record_list[[idx[1]]])
  } else {
    # Default to first record if no type strain available
    return(record_list[[1]])
  }
}

# Extract specific morphological traits from bacterial records
# Handles both single entries and lists of morphology blocks
extract_morphology_trait <- function(record, trait_name) {
  # Navigate to cell morphology section
  morph_block <- record$Morphology$`cell morphology`
  
  if (is.null(morph_block)) return("not reported")
  
  # Handle single-entry (not a list of blocks)
  # Direct access to trait if structure is flat
  if (!is.null(morph_block[[trait_name]])) {
    val <- morph_block[[trait_name]]
    if (is.null(val) || trimws(val) == "") return("not reported")
    return(trimws(val))
  }
  
  # Handle multi-entry list of blocks
  # Some records have multiple morphology observations
  if (is.list(morph_block) && all(sapply(morph_block, is.list))) {
    # Extract trait from each morphology block
    values <- unlist(lapply(morph_block, function(entry) {
      if (trait_name %in% names(entry)) entry[[trait_name]] else NULL
    }))
    
    # Clean up whitespace and remove empty values
    values <- trimws(values)
    values <- values[values != ""]
    
    if (length(values) == 0) {
      return("not reported")
    } else {
      # Combine unique values with semicolon separator
      return(paste(unique(values), collapse = "; "))
    }
  }
  
  return("not reported")
}

# Extract oxygen tolerance information from bacterial physiology data
# Oxygen tolerance indicates if bacteria are aerobic, anaerobic, etc.
extract_oxygen_tolerance <- function(record) {
  # Navigate to oxygen tolerance section in physiology data
  oxy_block <- record$`Physiology and metabolism`$`oxygen tolerance`
  
  if (is.null(oxy_block)) return("not reported")
  
  # Handle single-entry structure
  # Direct access to oxygen tolerance field
  if (!is.null(oxy_block$`oxygen tolerance`)) {
    values <- tolower(trimws(oxy_block$`oxygen tolerance`))
  } else if (is.list(oxy_block) && all(sapply(oxy_block, is.list))) {
    # Handle multiple oxygen tolerance entries
    # Extract from each block in the list
    values <- unlist(lapply(oxy_block, function(entry) {
      if ("oxygen tolerance" %in% names(entry)) {
        tolower(entry$`oxygen tolerance`)
      } else {
        NULL
      }
    }))
  } else {
    return("not reported")
  }
  
  # Remove duplicates and empty values
  values <- unique(values[values != ""])
  
  if (length(values) == 0) {
    return("not reported")
  } else {
    # Combine multiple values with semicolon separator
    return(paste(values, collapse = "; "))
  }
}

# Extract metabolites that bacteria can utilize (positive utilization)
# Returns vector of metabolite names where utilization activity is positive
extract_positive_metabolites <- function(record) {
  # Navigate to metabolite utilization data
  util <- record$`Physiology and metabolism`$`metabolite utilization`
  
  if (is.null(util) || length(util) == 0) return(character(0))
  
  # Filter for positive utilization activities only
  pos_metabolites <- unlist(lapply(util, function(entry) {
    # Check if entry has required fields and positive activity
    if (is.list(entry) &&
        "utilization activity" %in% names(entry) &&
        entry$`utilization activity` == "+" &&  # "+" indicates positive utilization
        "metabolite" %in% names(entry)) {
      return(entry$metabolite)
    } else {
      return(NULL)
    }
  }))
  
  # Return unique metabolite names, trimmed of whitespace
  unique(trimws(pos_metabolites))
}

# Extract enzymes with positive activity from bacterial records
# Returns vector of enzyme names where activity is positive
extract_positive_enzymes <- function(record) {
  # Navigate to enzyme data in physiology section
  enzymes <- record$`Physiology and metabolism`$enzymes
  
  if (is.null(enzymes) || length(enzymes) == 0) return(character(0))
  
  # Filter for positive enzyme activities only
  pos_enzymes <- unlist(lapply(enzymes, function(entry) {
    # Check if entry has required fields and positive activity
    if (is.list(entry) &&
        "activity" %in% names(entry) &&
        entry$activity == "+" &&  # "+" indicates positive enzyme activity
        "value" %in% names(entry)) {  # "value" contains the enzyme name
      return(entry$value)
    } else {
      return(NULL)
    }
  }))
  
  # Return unique enzyme names, trimmed of whitespace
  unique(trimws(pos_enzymes))
}


# Helper function to format traits
append_trait_rows <- function(species, set_type, set_ids, value = "+") {
  if (length(set_ids) == 0) return(data.frame())
  
  data.frame(
    species = species,
    set_type = set_type,
    set_id = set_ids,
    value = value,
    stringsAsFactors = FALSE
  )
}
