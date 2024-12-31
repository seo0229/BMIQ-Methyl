library(readxl)
library(ChAMP)
library(dplyr)
library(wateRmelon)

# Set file paths
methyl <- "/Users/elliottseo/Documents/GitHub/BMIQ-Methyl/Methyl_data/Mu EPIC KO v Wt AVG_beta & calculations-12.8.24.xlsx"
annot <- "/Users/elliottseo/Documents/GitHub/BMIQ-Methyl/Methyl_data/epic_mouse_annot_short_common_probes_manifests.xlsx"

# Loading methylation data
methyl_data <- read_excel(methyl, sheet = "Mu EPIC AVG_beta all CGs")

# Loading probe annotation data
annot_data <- read_excel(annot, sheet = "epic_mouse_annot_short_common_p")

# Extract probe type from targetid
annot_data <- annot_data %>%
    mutate(
        Probe_type = ifelse(substr(targetid, nchar(targetid) - 1, nchar(targetid) - 1) == "1", "Type1", "Type2")
    )

# Merge Methylation data with probe types
methyl_data <- methyl_data %>%
    rename(Probe = 'probe set') %>%
    left_join(dplyr::select(annot_data, targetid, Probe_type), by = c("Probe" = "targetid"))

# Ensure that the data is properly formatted for BMIQ
# - The beta values should be in a matrix format
# - Include only necessary columns for normalization
beta_values <- methyl_data %>%
  dplyr::select(-Probe, -Probe_type) %>%
  as.matrix()

# Check the distribution of probe types
probe_type_counts <- table(annot_data$Probe_type)
print(probe_type_counts)

# Ensure there are enough probes of each type
if (any(probe_type_counts < 2)) {
  stop("Not enough probes of each type for BMIQ normalization.")
}

# Create design vector for BMIQ
design.v <- ifelse(methyl_data$Probe_type == "Type1", 1, 2)

# Ensure design.v matches the rows of beta_values
design.v <- design.v[match(rownames(beta_values), methyl_data$Probe)]

# Run BMIQ normalization for each sample
normalized_data <- apply(beta_values, 2, function(beta) {
  wateRmelon::BMIQ(beta, design.v)$nbeta
})

# Convert the list of normalized data to a matrix
normalized_data <- do.call(cbind, normalized_data)

# Save normalized data to a CSV file
write.csv(normalized_data, "normalized_methylation_data.csv", row.names = TRUE)

# Optional: Snapshot the environment using renv
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)
renv::init()
renv::snapshot()