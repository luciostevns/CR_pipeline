#Libraries
library("tidyverse")


############### Initial filtering and wrangling ###############################

# Read the full IEDB dataset
autoimmune_data <- read_tsv("./Data/tcell_table_export_1740575887.tsv")

# Tidying data and filtering with relevant filters
autoimmune_data_wrangled <- autoimmune_data |>
  dplyr::rename(Assay_ID = "Assay ID - IEDB IRI",
                Sequence = `Epitope - Name`,
                Protein_source = `Epitope - Molecule Parent`,
                Protein_ID = `Epitope - Molecule Parent IRI`,
                Disease = `1st in vivo Process - Disease`,
                Disease_stage = `1st in vivo Process - Disease Stage`,
                MHC_restriction = `MHC Restriction - Name`,
                epitope_start_pos = `Epitope - Starting Position`,
                epitope_end_pos = `Epitope - Ending Position`) |>
  mutate(Sequence = str_trim(str_replace(Sequence, "\\+.*", "")),
         Protein_ID = str_extract(Protein_ID, "[^/]+$")
  ) |>
  filter(
    is.na(`Epitope - Modified residues`),
    str_length(Sequence) > 11 & str_length(Sequence) < 26
  ) |>
  distinct(Sequence, .keep_all = TRUE) |>
  dplyr::select(Assay_ID,Sequence, Protein_ID, Protein_source, Disease, Disease_stage, MHC_restriction,
                epitope_start_pos,epitope_end_pos)


########## Nested proteins removal #####################################


# Function to generate all possible 9-mers from a sequence
generate_9mers <- function(seq) {
  n <- nchar(seq)
  if (n < 9) return(character(0))  # Ignore sequences shorter than 9
  sapply(1:(n - 8), function(i) substr(seq, i, i + 8))
}

# Create a new column storing all 9-mers for each sequence
autoimmune_data_wrangled <- autoimmune_data_wrangled |>
  mutate(Nine_mers = lapply(Sequence, generate_9mers))

# Expand the dataframe so each row represents a 9-mer from the sequences
expanded_9mers <- autoimmune_data_wrangled |>
  dplyr::select(Sequence, Nine_mers) |>
  tidyr::unnest(Nine_mers)

# Identify sequences that are fully contained within longer ones
nested_sequences <- autoimmune_data_wrangled |>
  rowwise() |>
  mutate(Is_nested = any(expanded_9mers$Nine_mers %in% Nine_mers & 
                           nchar(expanded_9mers$Sequence) > nchar(Sequence))) |>
  ungroup()

# Keep only non-nested sequences
filtered_sequences <- nested_sequences |>
  filter(!Is_nested) |>
  dplyr::select(-Nine_mers, -Is_nested)  # Remove extra columns

write.csv(filtered_sequences, "./Data/wrangled_IEDB.csv", row.names = FALSE)

########################################################################################
