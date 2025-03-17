library(tidyverse)

temp_file <- tempfile(fileext = ".xlsx")
download.file("https://ictv.global/vmr/current", destfile = temp_file, mode = "wb")
vmr <- readxl::read_excel(temp_file, sheet = 2, col_types = "text")

# See all genome types
# vmr |>
#  select(Genome) |>
#  pull() |>
#  unique()

genome_types <- vmr |>
  select(Realm, Subrealm, Kingdom, Subkingdom, Phylum, Subphylum, Class, Subclass, Order, Suborder, Family, Subfamily, Genus, Subgenus, Species, Genome, `Exemplar or additional isolate`) |>
  filter(`Exemplar or additional isolate` == "E") |>
  mutate(Genome = case_when(
    Genome == "ssDNA(+/-)" ~ "ssDNA",
    Genome == "ssDNA(+)" ~ "ssDNA",
    Genome == "ssRNA(+/-)" ~ "ssRNA(-)",
    Genome == "ssRNA(-); ssRNA(+/-)" ~ "ssRNA(-)",
    Genome == "dsDNA-RT" ~ "dsDNA",
    Genome == "ssRNA-RT" ~ "ssRNA",
    Genome == "ssDNA(-)" ~ "ssDNA",
    T ~ Genome
  )) |>
  distinct()

vmr |>
  filter(Realm == "Riboviria" & Genome == "dsDNA")

# Define taxonomic levels in hierarchical order (from most to least inclusive)
tax_levels <- c(
  "Realm", "Subrealm", "Kingdom", "Subkingdom",
  "Phylum", "Subphylum", "Class", "Subclass",
  "Order", "Suborder", "Family", "Subfamily",
  "Genus", "Subgenus", "Species"
)

# Summarize the genome types for the realms with multiple genome types
realms <- genome_types |>
  select(Realm, Genome) |>
  filter(!is.na(Realm)) |>
  distinct() |>
  group_by(Realm) |>
  filter(n() > 1) |>
  mutate(Type = ifelse(str_detect(Genome, "DNA"), "DNA", "RNA")) |>
  summarize(Type = ifelse(n_distinct(Type) > 1, "uncharacterized", first(Type)), , Level = "Realm", .groups = "drop") |>
  rename(Taxon = Realm, Genome = Type)

# We'll store our summarized rows here.
summaries <- tibble()

# Copy the data so we can iteratively remove taxa already summarized.
remaining <- genome_types

# Loop over taxonomic levels:
for (lvl in tax_levels) {
  # For the current level, first keep only rows where that level is not NA.
  # Then group by that taxon and check whether all rows have the same Genome type.
  uniform_taxa <- remaining |>
    filter(!is.na(.data[[lvl]])) |>
    group_by(across(all_of(lvl))) |>
    # Only keep groups where Genome is uniform:
    filter(n_distinct(Genome) == 1) |>
    summarise(Genome = first(Genome), .groups = "drop") |>
    # Rename the grouping column as "Taxon" for clarity and record the taxonomic level.
    rename(Taxon = !!sym(lvl)) |>
    mutate(Level = lvl)

  # Append these uniform groups to our summaries.
  summaries <- bind_rows(summaries, uniform_taxa)

  # Remove all rows that belong to a taxon that has been summarized at this level.
  # (This ensures we don’t also output the lower-level rows from a group that is uniform.)
  remaining <- remaining |> filter(!(.data[[lvl]] %in% uniform_taxa$Taxon))
}

# Optionally, you can add any remaining (non-uniform) rows.
# (For example, if some viruses never reached a uniform group at any level.)
final_result <- bind_rows(realms, summaries, remaining)

# Write to TSV
final_result |>
  select(Taxon, Genome) |>
  rename(taxon = Taxon, pred_genome_type = Genome) |>
  write_tsv("suvtk/data/genome_types.tsv")
