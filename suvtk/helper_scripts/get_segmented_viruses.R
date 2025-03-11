library(tidyverse)

temp_file <- tempfile(fileext = ".xlsx")
download.file("https://ictv.global/vmr/current", destfile = temp_file, mode = "wb")
msl <- readxl::read_excel(temp_file, sheet = 2, col_types = "text")

# Calculate total counts for each taxonomic level from all records with a non-missing accession
kingdom_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Kingdom) %>%
  summarise(total = n(), .groups = "drop")

phylum_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Phylum) %>%
  summarise(total = n(), .groups = "drop")

subphylum_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Subphylum) %>%
  summarise(total = n(), .groups = "drop")

class_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Class) %>%
  summarise(total = n(), .groups = "drop")

order_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Order) %>%
  summarise(total = n(), .groups = "drop")

suborder_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Suborder) %>%
  summarise(total = n(), .groups = "drop")

family_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Family) %>%
  summarise(total = n(), .groups = "drop")

subfamily_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Subfamily) %>%
  summarise(total = n(), .groups = "drop")

genus_totals <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  group_by(Genus) %>%
  summarise(total = n(), .groups = "drop")

# Filter for segmented viruses: records with a colon and more than 1 segment
segmented_data <- msl %>%
  filter(!is.na(`Virus GENBANK accession`)) %>%
  filter(str_detect(`Virus GENBANK accession`, ":")) %>%
  mutate(segments = str_count(`Virus GENBANK accession`, ":"))


calculate_segment_stats <- function(df, level) {
  # Convert the string 'level' to a symbol for tidy evaluation
  level_sym <- sym(level)

  # Select the appropriate totals table based on the level provided
  totals <- switch(level,
    "Kingdom" = kingdom_totals,
    "Phylum" = phylum_totals,
    "Subphylum" = subphylum_totals,
    "Class" = class_totals,
    "Order" = order_totals,
    "Suborder" = suborder_totals,
    "Family" = family_totals,
    "Subfamily" = subfamily_totals,
    "Genus" = genus_totals,
    stop("Invalid level provided. Choose 'Order', 'Family', or 'Genus'.")
  )

  df %>%
    filter(segments > 1) %>%
    group_by(!!level_sym) %>%
    mutate(
      segmented = n(),
      min_segment = min(segments),
      max_segment = max(segments),
      majority_segment = as.integer(names(sort(table(segments), decreasing = TRUE)[1]))
    ) %>%
    ungroup() %>%
    select(!!level_sym, segmented, min_segment, max_segment, majority_segment) %>%
    # Join with the corresponding totals table. This assumes the totals table has a column with the same name.
    left_join(totals, by = setNames(level, level)) %>%
    mutate(segmented_fraction = segmented / total * 100) %>%
    rename(taxon = !!level_sym) %>%
    mutate(rank = level) %>%
    distinct()
}

# Kingdom-level segmented stats
kingdom_segmented <- calculate_segment_stats(segmented_data, "Kingdom")

# Phylum-level segmented stats
phylum_segmented <- calculate_segment_stats(segmented_data, "Phylum")

# Phylum-level segmented stats
subphylum_segmented <- calculate_segment_stats(segmented_data, "Subphylum")

# Class-level segmented stats
class_segmented <- calculate_segment_stats(segmented_data, "Class")

# Order-level segmented stats
order_segmented <- calculate_segment_stats(segmented_data, "Order")

# Suborder-level segmented stats
suborder_segmented <- calculate_segment_stats(segmented_data, "Suborder")

# Family-level segmented stats
family_segmented <- calculate_segment_stats(segmented_data, "Family")

# Subfamily-level segmented stats
subfamily_segmented <- calculate_segment_stats(segmented_data, "Subfamily")

# Genus-level segmented stats
genus_segmented <- calculate_segment_stats(segmented_data, "Genus")

# Create mapping tables for parent relationships
# For Class: use Subphylum if available, otherwise Phylum
class_parent_map <- msl %>%
  select(Class, Subphylum, Phylum) %>%
  distinct()

# For Family: use Suborder if available, otherwise Order
family_parent_map <- msl %>%
  select(Family, Suborder, Order) %>%
  distinct()

# For Genus: use Subfamily if available, otherwise Family
genus_parent_map <- msl %>%
  select(Genus, Subfamily, Family) %>%
  distinct()

# Add parent column for Kingdom, Phylum and Order
kingdom_map <- msl %>%
  select(Kingdom, Realm) %>%
  distinct()

phylum_map <- msl %>%
  select(Phylum, Kingdom) %>%
  distinct()

order_map <- msl %>%
  select(Order, Class) %>%
  distinct()

# For Subphylum, Suborder and Subfamily, we'll assign their parent as follows:
# For Subphylum: parent = Phylum
# For Suborder: parent = Order
# For Subfamily: parent = Family
subphylum_map <- msl %>%
  select(Subphylum, Phylum) %>%
  distinct()

suborder_map <- msl %>%
  select(Suborder, Order) %>%
  distinct()

subfamily_map <- msl %>%
  select(Subfamily, Family) %>%
  distinct()

# Add parent column for each level
kingdom_segmented <- kingdom_segmented %>%
  left_join(kingdom_map, by = c("taxon" = "Kingdom")) %>%
  rename(parent = Realm)

phylum_segmented <- phylum_segmented %>%
  left_join(phylum_map, by = c("taxon" = "Phylum")) %>%
  rename(parent = Kingdom)

subphylum_segmented <- subphylum_segmented %>%
  left_join(subphylum_map, by = c("taxon" = "Subphylum")) %>%
  rename(parent = Phylum)

# For Class: Use Subphylum if available; if NA, then use Phylum.
class_segmented <- class_segmented %>%
  left_join(class_parent_map, by = c("taxon" = "Class")) %>%
  mutate(parent = if_else(!is.na(Subphylum) & Subphylum != "", Subphylum, Phylum)) %>%
  select(-Subphylum, -Phylum)

order_segmented <- order_segmented %>%
  left_join(order_map, by = c("taxon" = "Order")) %>%
  rename(parent = Class)

suborder_segmented <- suborder_segmented %>%
  left_join(suborder_map, by = c("taxon" = "Suborder")) %>%
  rename(parent = Order)

# For Family: Use Suborder if available; if NA, then use Order.
family_segmented <- family_segmented %>%
  left_join(family_parent_map, by = c("taxon" = "Family")) %>%
  mutate(parent = if_else(!is.na(Suborder) & Suborder != "", Suborder, Order)) %>%
  select(-Suborder, -Order)

# For Subfamily: Parent is Family (as typical)
subfamily_segmented <- subfamily_segmented %>%
  left_join(subfamily_map, by = c("taxon" = "Subfamily")) %>%
  rename(parent = Family)

# For Genus: Use Subfamily if available; if NA, then use Family.
genus_segmented <- genus_segmented %>%
  left_join(genus_parent_map, by = c("taxon" = "Genus")) %>%
  mutate(parent = if_else(!is.na(Subfamily) & Subfamily != "", Subfamily, Family)) %>%
  select(-Subfamily, -Family)

# Combine all levels into one long-format data frame
final_long <- bind_rows(
  kingdom_segmented, phylum_segmented, subphylum_segmented, class_segmented,
  order_segmented, suborder_segmented, family_segmented, subfamily_segmented, genus_segmented
) %>%
  select(
    rank, taxon, parent, total, segmented, segmented_fraction,
    majority_segment, min_segment, max_segment
  ) %>%
  filter(!is.na(taxon))

# Write to TSV
write_tsv(final_long, "gb_submitter/data/segmented_viruses.tsv")
