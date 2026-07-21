# Heatmaps of monthly PHI-base pathogen detections
# Written July 2026

# Using the filtered dfs created in the Host_Bar_Charts.R script

# Packages --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(taxize)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(patchwork)


# Settings --------------------------------------------------------------------

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")

# Remove phyla represented by fewer than this number of detected pathogen taxa
minimum_pathogens_per_phylum <- 2

# ColourBrewer sequential palette for heatmaps
heatmap_colours <- brewer.pal(
  n = 9,
  name = "YlGnBu"
)

# Heatmap theme 

heatmap_theme <- theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
  )

# Read processed data ---------------------------------------------------------

all_data_filtered <- readRDS(
  "Derived_Phibase_Data/all_data_filtered.rds"
)

plot_data <- readRDS(
  "Derived_Phibase_Data/plant_host_plot_data.rds"
)

# Monthly heatmap: Fragaria versus other plant hosts ----------------------------------
## Pathogen number detected ------------------------------------------------
host_category_month_counts <- plot_data %>%
  filter(!is.na(month), !is.na(location), !is.na(taxaID)
  ) %>%
  mutate(
    host_category = if_else(
      is_fragaria,
      "Fragaria-associated",
      "Other plant hosts"
    )
  ) %>%
  distinct(
    month,
    location,
    host_category,
    taxaID
  ) %>%
  count(
    month,
    location,
    host_category,
    name = "n_pathogens"
  ) %>%
  complete(
    month,
    location,
    host_category,
    fill = list(
      n_pathogens = 0
    )
  ) %>%
  mutate(
    host_category = factor(
      host_category,
      levels = c(
        "Fragaria-associated",
        "Other plant hosts"
      )
    )
  )

### Plot pathogen number ---------
# With a shared gradient, this makes fragaria look low as there are only 6 potential species 
host_category_heatmap <- ggplot(
  host_category_month_counts,
  aes(
    x = month,
    y = location,
    fill = n_pathogens
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ host_category,
    nrow = 1
  ) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    breaks = pretty_breaks(),
    name = "Pathogen\nspecies (n)"
  ) +
  labs(
    x = "Month",
    y = "Collection location"
  ) +
  theme_minimal(base_size = 10)

host_category_heatmap

### Plot separate and combine ------------------------
fragaria_month_counts <- host_category_month_counts %>%
  filter(host_category == "Fragaria-associated")

# Fragaria hosts
fragaria_heatmap <- ggplot(
  fragaria_month_counts,
  aes(
    x = month,
    y = location,
    fill = n_pathogens
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
        label = n_pathogens,
        colour = if_else(
          n_pathogens >= 4,
          "white",
          "black"
        )
    ),
    size = 3
  ) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    limits = c(0, 6),
    breaks = 0:6,
    name = "Pathogen\nspecies (n)"
  ) +
  scale_colour_identity() +
  labs(
    title = expression(italic("Fragaria") * "-associated"),
    x = "Month",
    y = "Collection location"
  ) +
  heatmap_theme


# Other hosts
other_host_month_counts <- host_category_month_counts %>%
  filter(host_category == "Other plant hosts")

other_host_heatmap <- ggplot(
  other_host_month_counts,
  aes(
    x = month,
    y = location,
    fill = n_pathogens
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    breaks = pretty_breaks(),
    name = "Pathogen\nspecies (n)"
  ) +
  labs(
    title = "Other plant hosts",
    x = "Month",
    y = NULL
  ) +
  heatmap_theme

### Combine and export -------------
host_category_heatmap <- fragaria_heatmap + other_host_heatmap +
  plot_layout(
    widths = c(1, 1),
    guides = "keep" # Each have seperate scale
  )

host_category_heatmap

ggsave(
  filename = "Graphs/PHIbase_graphs/PHIbase_Heatmap_Fragaria_Other.pdf",
  plot = host_category_heatmap,
  width = 172,
  height = 100,
  units = "mm",
  device = cairo_pdf
)


## Normalised pathogen abundance by host category -----------------------------

# A pathogen associated with both Fragaria and another plant host will
# contribute to both categories.
#
# distinct() prevents multiple PHI-base host records from duplicating the
# pathogen's normalised read count within the same category.

sample_pathogen_categories <- plot_data %>%
  mutate(
    host_category = if_else(
      is_fragaria,
      "Fragaria-associated",
      "Other plant hosts"
    )
  ) %>%
  distinct(
    Sample,
    month,
    location,
    taxaID,
    host_category,
    .keep_all = TRUE
  )

sample_host_abundance <- sample_pathogen_categories %>%
  group_by(
    Sample,
    month,
    location,
    host_category
  ) %>%
  summarise(
    normalised_hits = sum(
      hits_per_100000,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# At present, a sample only appears in a category when at least one relevant pathogen was detected. 
# To calculate a proper mean, samples with no detections in that category need to contribute zero.
sample_metadata <- all_data_filtered %>%
  distinct(
    Sample,
    month,
    location
  )

sample_host_abundance_complete <- sample_metadata %>%
  tidyr::crossing(
    host_category = c(
      "Fragaria-associated",
      "Other plant hosts"
    )
  ) %>%
  left_join(
    sample_host_abundance,
    by = c(
      "Sample",
      "month",
      "location",
      "host_category"
    )
  ) %>%
  mutate(
    normalised_hits = replace_na(
      normalised_hits,
      0
    ),
    host_category = factor(
      host_category,
      levels = c(
        "Fragaria-associated",
        "Other plant hosts"
      )
    )
  )

# Average normalised abundance in each month and location
host_category_month_abundance <- sample_host_abundance_complete %>%
  group_by(
    month,
    location,
    host_category
  ) %>%
  summarise(
    mean_normalised_hits = mean(
      normalised_hits,
      na.rm = TRUE
    ),
    median_normalised_hits = median(
      normalised_hits,
      na.rm = TRUE
    ),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  )

# Separate for the heatmap - Not sure I like logged or un-logged of this
fragaria_month_abundance <- host_category_month_abundance %>%
  filter(
    host_category == "Fragaria-associated"
  )

other_host_month_abundance <- host_category_month_abundance %>%
  filter(
    host_category == "Other plant hosts"
  )

fragaria_abundance_heatmap <- ggplot(
  fragaria_month_abundance,
  aes(
    x = month,
    y = location,
    fill = mean_normalised_hits
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  # Logged
  scale_fill_gradientn(
    colours = heatmap_colours,
    trans = scales::pseudo_log_trans(
      sigma = 1
    ),
    name = "Mean hits per\n100,000 reads"
  ) +
  # Not transformed
  # scale_fill_distiller(
  #   palette = "YlGnBu",
  #   direction = 1,
  #   breaks = c(
  #     0,
  #     5000,
  #     10000,
  #     15000,
  #     20000,
  #     25000
  #   ),
  #   labels = scales::label_number(
  #     big.mark = ","
  #   ),
  #   name = "Mean hits per\n100,000 reads"
  # ) +
  labs(
    title = expression(
      italic("Fragaria") * "-associated"
    ),
    x = "Month",
    y = "Collection location"
  ) +
  heatmap_theme

fragaria_abundance_heatmap

# other host heatmap
other_host_abundance_heatmap <- ggplot(
  other_host_month_abundance,
  aes(
    x = month,
    y = location,
    fill = mean_normalised_hits
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    trans = scales::pseudo_log_trans(
      sigma = 1
    ),
    name = "Mean hits per\n100,000 reads"
  ) +
  labs(
    title = "Other plant hosts",
    x = "Month",
    y = NULL
  ) +
  heatmap_theme

other_host_abundance_heatmap

### Combine and export --------
host_abundance_heatmap <- (
  fragaria_abundance_heatmap +
    other_host_abundance_heatmap
) +
  patchwork::plot_layout(
    widths = c(1, 1),
    guides = "keep"
  )

host_abundance_heatmap

ggsave(
  filename = "Graphs/PHIbase_graphs/PHIbase_Heatmap_norm_Fragaria_Other.pdf",
  plot = host_abundance_heatmap,
  width = 172,
  height = 100,
  units = "mm",
  device = cairo_pdf
)


# Phylum level pathogen heatmap -----------------------------------------------

## Pathogen taxonomy lookup ----------------------------------------------------

pathogen_taxonomy_file <- paste0(
  "Derived_Phibase_Data/",
  "pathogen_taxonomy_lookup.rds"
)


# Function to extract one taxonomic rank

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
    filter(
      rank == target_rank
    ) %>%
    pull(name) %>%
    first(
      default = NA_character_
    )
}


# IDs currently present in the filtered sequencing data

pathogen_ids <- all_data_filtered %>%
  distinct(taxaID) %>%
  arrange(taxaID)


# Load the existing taxonomy cache where available

if (file.exists(pathogen_taxonomy_file)) {
  
  pathogen_taxonomy_df <- readRDS(
    pathogen_taxonomy_file
  ) %>%
    mutate(
      taxaID = as.character(taxaID)
    )
  
} else {
  
  pathogen_taxonomy_df <- tibble(
    taxaID = character(),
    pathogen_phylum = character()
  )
}

# Identify IDs that have not yet been queried
missing_pathogen_ids <- pathogen_ids %>%
  anti_join(
    pathogen_taxonomy_df,
    by = "taxaID"
  )

# Query only previously unseen pathogen IDs

if (nrow(missing_pathogen_ids) > 0) {
  
  pathogen_ids_to_query <- missing_pathogen_ids %>%
    pull(taxaID)
  
  pathogen_taxonomy <- taxize::classification(
    pathogen_ids_to_query,
    db = "ncbi"
  )
  
  taxonomy_list <- pathogen_taxonomy %>%
    unclass() %>%
    as.list() %>%
    unname()
  
  
  if (
    length(taxonomy_list) !=
    length(pathogen_ids_to_query)
  ) {
    stop(
      paste(
        "The number of pathogen taxonomy results does not",
        "match the number of taxonomy IDs queried."
      )
    )
  }
  
  
  new_pathogen_taxonomy <- tibble(
    taxaID = pathogen_ids_to_query,
    taxonomy_result = taxonomy_list
  ) %>%
    mutate(
      pathogen_phylum = map_chr(
        taxonomy_result,
        extract_taxonomic_rank,
        target_rank = "phylum"
      )
    ) %>%
    select(
      taxaID,
      pathogen_phylum
    )
  
  pathogen_taxonomy_df <- pathogen_taxonomy_df %>%
    bind_rows(
      new_pathogen_taxonomy
    ) %>%
    distinct(
      taxaID,
      .keep_all = TRUE
    )
  
  saveRDS(
    pathogen_taxonomy_df,
    pathogen_taxonomy_file
  )
}


# Inspect the taxonomy lookup

pathogen_taxonomy_df %>%
  count(
    pathogen_phylum,
    sort = TRUE
  )

## Add pathogen phylum to sequencing data --------------------------------------

all_data_with_pathogen_taxonomy <- all_data_filtered %>%
  left_join(
    pathogen_taxonomy_df,
    by = "taxaID"
  ) %>%
  mutate(
    pathogen_phylum = replace_na(
      pathogen_phylum,
      "Phylum unavailable"
    )
  )

# Determine how many pathogen taxa belong to each phylum ----------------------

phylum_totals <- all_data_with_pathogen_taxonomy %>%
  distinct(
    pathogen_phylum,
    taxaID
  ) %>%
  count(
    pathogen_phylum,
    name = "n_pathogens_total"
  ) %>%
  arrange(
    desc(n_pathogens_total)
  )


# Keep phyla represented by at least the selected number of pathogen taxa
phyla_to_plot <- phylum_totals %>%
  filter(
    n_pathogens_total >= minimum_pathogens_per_phylum
  )


# Inspect phyla excluded from the heatmap (4)
phyla_removed <- phylum_totals %>%
  filter(
    n_pathogens_total < minimum_pathogens_per_phylum
  )
phyla_removed

## Count unique pathogen taxa by month, location and phylum --------------------

phylum_month_counts <- all_data_with_pathogen_taxonomy %>%
  filter(
    !is.na(month),
    !is.na(location)
  ) %>%
  semi_join(
    phyla_to_plot,
    by = "pathogen_phylum"
  ) %>%
  distinct(
    month,
    location,
    pathogen_phylum,
    taxaID
  ) %>%
  count(
    month,
    location,
    pathogen_phylum,
    name = "n_pathogens"
  ) %>%
  complete(
    month,
    location,
    pathogen_phylum,
    fill = list(
      n_pathogens = 0
    )
  )


# Order phyla by their overall number of detected taxa
phylum_order <- phyla_to_plot %>%
  arrange(
    desc(n_pathogens_total)
  ) %>%
  pull(pathogen_phylum)


phylum_month_counts <- phylum_month_counts %>%
  mutate(
    pathogen_phylum = factor(
      pathogen_phylum,
      levels = rev(phylum_order)
    )
  )

## Plot phylum heatmap counts
phylum_heatmap <- ggplot(
phylum_month_counts,
aes(
  x = month,
  y = pathogen_phylum,
  fill = n_pathogens
)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  facet_wrap(
    ~ location,
    nrow = 1
  ) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    breaks = pretty_breaks(),
    name = "Pathogen\nspecies (n)"
  ) +
  labs(
    x = "Month",
    y = "Pathogen phylum"
  ) +
  heatmap_theme

phylum_heatmap

## Export phylum heatmap counts
ggsave(
  filename = "Graphs/PHIbase_graphs/PHIbase_Heatmap_Phylum.pdf",
  plot = phylum_heatmap,
  width = 172,
  height = 100,
  units = "mm",
  device = cairo_pdf
)
