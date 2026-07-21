# Script to get host taxonomy for PHI-base minimap alignment data
# Written July 2026


# Packages --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(taxize)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Settings --------------------------------------------------------------------

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")

# Families represented by fewer pathogen species than this will not be plotted
minimum_pathogens_per_family <- 2

plant_host_categories <- c(
  "eudicots",
  "monocots",
  "seed plants",
  "flowering plants"
)

#  Plot theme
custom_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

# Colours
location_colours <- setNames(
  brewer.pal(3, "Dark2"),
  c("Field", "Greenhouse 1", "Greenhouse 2")
)

# Alignment data --------------------------------------------------------------

# Read the main alignment dataset

main_alignment_data <- read.csv(
  "Tiptree_taxa_counts_initial_analysis/Q20_Tiptree_combined_data.csv",
  header = TRUE
) %>%
  select(
    taxaID,
    read_count,
    species,
    barcode,
    total_reads,
    Sample,
    Date.collected,
    repeat.,
    location
  )


# Read the December alignment dataset

december_alignment_data <- read.csv(
  "Tiptree_taxa_counts_initial_analysis/Tiptree_Dec_unfiltered.csv",
  header = TRUE
) %>%
  select(
    all_of(names(main_alignment_data))
  )


# Combine datasets

all_data <- bind_rows(
  main_alignment_data,
  december_alignment_data
)


## Process sample information --------------------------------------------------

all_data <- all_data %>%
  mutate(
    taxaID = as.character(taxaID),
    
    Date.collected = as.Date(
      as.character(Date.collected),
      format = "%d/%m/%Y"
    ),
    
    month = factor(
      format(Date.collected, "%b"),
      levels = c(
        "Dec",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug"
      )
    ),
    
    location = factor(
      as.character(location),
      levels = c("F", "RG", "LG"),
      labels = c(
        "Field",
        "Greenhouse 1",
        "Greenhouse 2"
      )
    ),
    
    Sample = as.character(Sample)
  )


## Filter and normalise sequencing data ----------------------------------------

all_data_filtered <- all_data %>%
  filter(
    read_count >= 10,
    !is.na(taxaID),
    taxaID != "",
    !is.na(total_reads),
    total_reads > 0
  ) %>%
  mutate(
    hits_per_100000 = read_count * 100000 / total_reads,
    log_hits = log10(hits_per_100000)
  )


# Check the filtered alignment data
alignment_summary <- all_data_filtered %>%
  summarise(
    n_rows = n(),
    n_samples = n_distinct(Sample),
    n_pathogen_taxa = n_distinct(taxaID),
    n_locations = n_distinct(location, na.rm = TRUE)
  )

alignment_summary


# PHI-base data ---------------------------------------------------------------

phibase <- read.csv(
  paste0(
    "~/OneDrive - Norwich BioScience Institutes/",
    "Pathogen_Database/Pathogen_Database_111224/",
    "phibase_4.17.csv"
  )
)


# Optional check of host categories in PHI-base (This is how I selected those to filter at the start)

# phibase_host_category_counts <- phibase %>%
#   count(
#     Host_descripton,
#     sort = TRUE
#   )

## Detected pathogen IDs -------------------------------------------------------

detected_pathogen_ids <- all_data_filtered %>%
  distinct(taxaID)


## Create the PHI-base plant-host lookup ---------------------------------------

# Each row represents one unique pathogen-host association.
# A pathogen may therefore appear more than once if it has multiple hosts.

phibase_plant_hosts <- phibase %>%
  transmute(
    taxaID = as.character(
      Pathogen_NCBI_species_Taxonomy.ID
    ),
    
    host_taxaID = as.character(
      Host_NCBI_Taxonomy_ID
    ),
    
    host_description = Host_descripton
  ) %>%
  filter(
    !is.na(taxaID),
    taxaID != "",
    !is.na(host_taxaID),
    host_taxaID != "",
    host_description %in% plant_host_categories
  ) %>%
  distinct(
    taxaID,
    host_taxaID,
    .keep_all = TRUE
  ) %>%
  semi_join(
    detected_pathogen_ids,
    by = "taxaID"
  ) %>%
  mutate(
    host_source = "PHI-base"
  )


## Manually added pathogen-host associations -----------------------------------
# These aren't in PHIbase with correct host

manual_pathogen_hosts <- tibble(
  taxaID = c(
    "79252",  # Podosphaera aphanis
    "5025"    # Venturia inaequalis
  ),
  
  host_taxaID = c(
    "3747",   # Fragaria x ananassa
    "3750"    # Malus domestica
  ),
  
  host_description = c(
    "eudicots",
    "eudicots"
  ),
  
  host_source = "Manually added"
)


# Check whether the manually added pathogens occur in the alignment data

manual_pathogens_not_detected <- manual_pathogen_hosts %>%
  anti_join(
    detected_pathogen_ids,
    by = "taxaID"
  )

manual_pathogens_not_detected


## Combine PHI-base and manually added associations ----------------------------

phibase_hosts_complete <- phibase_plant_hosts %>%
  bind_rows(
    manual_pathogen_hosts
  ) %>%
  semi_join(
    detected_pathogen_ids,
    by = "taxaID"
  ) %>%
  distinct(
    taxaID,
    host_taxaID,
    .keep_all = TRUE
  )


# Check for duplicated pathogen-host pairs

duplicate_host_pairs <- phibase_hosts_complete %>%
  count(
    taxaID,
    host_taxaID
  ) %>%
  filter(n > 1)

duplicate_host_pairs


# Summarise plant-host associations before taxonomy lookup

plant_host_lookup_summary <- phibase_hosts_complete %>%
  summarise(
    n_pathogens_with_plant_hosts = n_distinct(taxaID),
    n_unique_host_taxa = n_distinct(host_taxaID),
    n_pathogen_host_pairs = n()
  )

plant_host_lookup_summary


## Record detected pathogens excluded from the plant-host analysis -------------

# This is an audit table, not an error check.
# These pathogens were detected but do not have a retained plant host.

pathogens_without_plant_hosts <- detected_pathogen_ids %>%
  anti_join(
    phibase_hosts_complete %>%
      distinct(taxaID),
    by = "taxaID"
  ) %>%
  left_join(
    all_data_filtered %>%
      distinct(
        taxaID,
        species
      ),
    by = "taxaID"
  )

pathogens_without_plant_hosts


# Retrieve host taxonomy from NCBI --------------------------------------------

host_id_lookup <- phibase_hosts_complete %>%
  distinct(host_taxaID) %>%
  arrange(host_taxaID)

host_ids <- host_id_lookup %>%
  pull(host_taxaID)


host_taxonomy <- taxize::classification(
  host_ids,
  db = "ncbi"
)


# Convert the special classification object into an ordinary list

taxonomy_list <- host_taxonomy %>%
  unclass() %>%
  as.list() %>%
  unname()


# Stop if the number of taxonomy results does not match the number queried

if (length(taxonomy_list) != length(host_ids)) {
  stop(
    paste(
      "The number of NCBI taxonomy results does not match",
      "the number of host taxonomy IDs queried."
    )
  )
}


# Pair each taxonomy result with its original host ID

taxonomy_lookup <- tibble(
  host_taxaID = host_ids,
  taxonomy_result = taxonomy_list
)


# Identify failed taxonomy queries

failed_host_ids <- taxonomy_lookup %>%
  filter(
    map_lgl(
      taxonomy_result,
      ~ is.null(.x) ||
        !is.data.frame(.x) ||
        nrow(.x) == 0
    )
  ) %>%
  select(host_taxaID)

failed_host_ids


## Function to extract one taxonomic rank --------------------------------------

extract_taxonomic_rank <- function(
    classification_table,
    target_rank
) {
  
  if (
    is.null(classification_table) ||
    !is.data.frame(classification_table) ||
    nrow(classification_table) == 0
  ) {
    return(NA_character_)
  }
  
  classification_table %>%
    filter(rank == target_rank) %>%
    pull(name) %>%
    first(default = NA_character_)
}


## Convert taxonomy results into a dataframe -----------------------------------

host_taxonomy_df <- taxonomy_lookup %>%
  mutate(
    host_species = map_chr(
      taxonomy_result,
      extract_taxonomic_rank,
      target_rank = "species"
    ),
    
    host_genus = map_chr(
      taxonomy_result,
      extract_taxonomic_rank,
      target_rank = "genus"
    ),
    
    host_family = map_chr(
      taxonomy_result,
      extract_taxonomic_rank,
      target_rank = "family"
    ),
    
    host_order = map_chr(
      taxonomy_result,
      extract_taxonomic_rank,
      target_rank = "order"
    )
  ) %>%
  select(
    host_taxaID,
    host_species,
    host_genus,
    host_family,
    host_order
  )


# Check the taxonomy extraction

taxonomy_summary <- host_taxonomy_df %>%
  summarise(
    n_hosts = n(),
    n_missing_species = sum(is.na(host_species)),
    n_missing_genus = sum(is.na(host_genus)),
    n_missing_family = sum(is.na(host_family)),
    n_missing_order = sum(is.na(host_order))
  )

taxonomy_summary


# Inspect hosts for which no family was retrieved

hosts_without_family <- host_taxonomy_df %>%
  filter(
    is.na(host_family)
  )

hosts_without_family


# Add taxonomy to the pathogen-host lookup ------------------------------------

phibase_hosts_taxonomy <- phibase_hosts_complete %>%
  left_join(
    host_taxonomy_df,
    by = "host_taxaID"
  ) %>%
  mutate(
    is_fragaria = replace_na(
      host_genus == "Fragaria",
      FALSE
    )
  )


# Inspect the completed pathogen-host lookup

phibase_hosts_taxonomy %>%
  arrange(
    taxaID,
    host_family,
    host_genus,
    host_species
  )


## Summarise the host range of each detected pathogen --------------------------

pathogen_host_summary <- phibase_hosts_taxonomy %>%
  group_by(taxaID) %>%
  summarise(
    n_host_taxa = n_distinct(host_taxaID),
    
    n_host_families = n_distinct(
      host_family,
      na.rm = TRUE
    ),
    
    host_families = paste(
      sort(unique(host_family[!is.na(host_family)])),
      collapse = "; "
    ),
    
    associated_with_fragaria = any(is_fragaria),
    
    .groups = "drop"
  )

pathogen_host_summary


# Join host information to sequencing data -----------------------------------

# inner_join() deliberately removes detected pathogens without a retained
# plant-host association.
#
# The join is many-to-many because:
# 1. each pathogen may be detected in several sequencing samples;
# 2. each pathogen may have several host associations.

plot_data <- all_data_filtered %>%
  inner_join(
    phibase_hosts_taxonomy,
    by = "taxaID",
    relationship = "many-to-many"
  ) %>%
  select(
    # Sampling information
    Sample,
    barcode,
    Date.collected,
    month,
    repeat.,
    location,
    
    # Pathogen and sequencing information
    taxaID,
    species,
    read_count,
    total_reads,
    hits_per_100000,
    log_hits,
    
    # Plant-host information
    host_taxaID,
    host_species,
    host_genus,
    host_family,
    host_order,
    host_description,
    host_source,
    is_fragaria
  ) %>%
  arrange(
    location,
    Date.collected,
    Sample,
    taxaID,
    host_family,
    host_species
  )


# Final check of the joined plotting data -------------------------------------

plot_data_summary <- plot_data %>%
  summarise(
    n_rows = n(),
    n_samples = n_distinct(Sample),
    n_pathogens = n_distinct(taxaID),
    n_host_taxa = n_distinct(host_taxaID),
    n_host_families = n_distinct(
      host_family,
      na.rm = TRUE
    ),
    n_fragaria_pathogens = n_distinct(
      taxaID[is_fragaria]
    ),
    n_locations = n_distinct(location)
  )

plot_data_summary


# Check for unexpected missing host information

plot_data_missing_values <- plot_data %>%
  summarise(
    n_missing_host_ids = sum(is.na(host_taxaID)),
    n_missing_host_families = sum(is.na(host_family)),
    n_non_plant_categories = sum(
      !host_description %in% plant_host_categories
    )
  )

plot_data_missing_values


# Family-level plotting data --------------------------------------------------

# Each pathogen is counted once per host family overall.
# A pathogen associated with more than one family is counted in each family.

all_family_totals <- plot_data %>%
  filter(
    !is.na(host_family)
  ) %>%
  distinct(
    host_family,
    taxaID
  ) %>%
  count(
    host_family,
    name = "n_pathogens_total"
  ) %>%
  arrange(
    desc(n_pathogens_total)
  )


# Families retained for the plot

family_totals <- all_family_totals %>%
  filter(
    n_pathogens_total >= minimum_pathogens_per_family
  )


# Families excluded from the plot

families_removed <- all_family_totals %>%
  filter(
    n_pathogens_total < minimum_pathogens_per_family
  )

families_removed


# Count unique detected pathogens in each family and location

family_location_counts <- plot_data %>%
  filter(
    !is.na(location),
    !is.na(host_family)
  ) %>%
  semi_join(
    family_totals,
    by = "host_family"
  ) %>%
  distinct(
    location,
    host_family,
    taxaID
  ) %>%
  count(
    host_family,
    location,
    name = "n_pathogens"
  ) %>%
  complete(
    host_family,
    location,
    fill = list(
      n_pathogens = 0
    )
  ) %>%
  mutate(
    host_family = factor(
      host_family,
      levels = family_totals$host_family
    )
  )


# Inspect the dataframe used for the family plot

family_location_counts


## Family-level plot -----------------------------------------------

# This is deliberately simple for now.
# We can make it publication-ready in the next step.

family_plot <- ggplot(
  family_location_counts,
  aes(
    x = host_family,
    y = n_pathogens,
    fill = location
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = location_colours,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = "Host family",
    y = "Pathogen species (n)",
    fill = "Collection location"
  ) +
  custom_theme

family_plot


# Fragaria-specific plotting data ---------------------------------------------

## Fragaria location plot  ---------------------------------------------

# Count each Fragaria-associated pathogen once at each location,
# regardless of how many samples it was detected in there.

fragaria_location_counts <- plot_data %>%
  filter(
    is_fragaria,
    !is.na(location),
    !is.na(taxaID)
  ) %>%
  distinct(
    location,
    taxaID
  ) %>%
  count(
    location,
    name = "n_fragaria_pathogens"
  ) %>%
  complete(
    location,
    fill = list(
      n_fragaria_pathogens = 0
    )
  )

fragaria_location_counts

fragaria_plot <- ggplot(
  fragaria_location_counts,
  aes(
    x = location,
    y = n_fragaria_pathogens,
    fill = location
  )
) +
  geom_col(
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = location_colours,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = "Collection location",
    y = expression(
      paste(italic("Fragaria"), "-associated pathogens (n)")
    )
  ) +
  custom_theme 

fragaria_plot

# Combined plot for publication ---------
combined_plot <- family_plot + fragaria_plot +
  plot_layout(
    widths = c(8, 1.5),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right"
  )

combined_plot

ggsave(
  filename = "~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree/Graphs/PHIbase_host_family_and_Fragaria.pdf",
  plot = combined_plot,
  width = 172,
  height = 100,
  units = "mm",
  device = cairo_pdf
)

# Unwanted plot ----------------
## Fragaria species plots -----------------------------------------
# A pathogen may be linked to more than one Fragaria host record.
# Keep only one copy of each pathogen in each sequencing sample.

fragaria_sample_data <- plot_data %>%
  filter(
    is_fragaria
  ) %>%
  distinct(
    Sample,
    location,
    taxaID,
    .keep_all = TRUE
  ) %>%
  mutate(
    pathogen_label = if_else(
      !is.na(species) & species != "",
      species,
      taxaID
    )
  )

# Check which Fragaria-associated pathogens were detected

fragaria_pathogens_detected <- fragaria_sample_data %>%
  distinct(
    taxaID,
    pathogen_label
  ) %>%
  arrange(pathogen_label)

fragaria_pathogens_detected

# Summarise Fragaria-associated detections by location 

# n_samples_detected records how many samples at each location contained
# each pathogen.
#
# The abundance summaries are also retained in case you later decide to plot
# normalised read counts rather than sample-level detection frequency.

fragaria_location_summary <- fragaria_sample_data %>%
  group_by(
    location,
    taxaID,
    pathogen_label
  ) %>%
  summarise(
    n_samples_detected = n_distinct(Sample),
    
    mean_hits_per_100000 = mean(
      hits_per_100000,
      na.rm = TRUE
    ),
    
    median_hits_per_100000 = median(
      hits_per_100000,
      na.rm = TRUE
    ),
    
    maximum_hits_per_100000 = max(
      hits_per_100000,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# Order Fragaria-associated pathogens by total detection frequency

fragaria_pathogen_order <- fragaria_location_summary %>%
  group_by(pathogen_label) %>%
  summarise(
    total_samples_detected = sum(n_samples_detected),
    .groups = "drop"
  ) %>%
  arrange(
    desc(total_samples_detected)
  ) %>%
  pull(pathogen_label)


fragaria_location_summary <- fragaria_location_summary %>%
  mutate(
    pathogen_label = factor(
      pathogen_label,
      levels = rev(fragaria_pathogen_order)
    )
  )


# Inspect the dataframe for the Fragaria plot

fragaria_location_summary


# Plot

# This currently shows the number of samples in which each pathogen was
# detected at each location.

fragaria_plot <- ggplot(
  fragaria_location_summary,
  aes(
    x = pathogen_label,
    y = n_samples_detected,
    fill = location
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(
    x = "Fragaria-associated pathogen",
    y = "Number of samples detected",
    fill = "Collection location"
  )  + custom_theme

fragaria_plot
