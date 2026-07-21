# Strawberry Analysis

This repository contains scripts used to analyse sequencing data generated from air sampling at a strawberry farm. The analysis combines pathogen sequencing outputs with sample metadata, DNA yield measurements, disease scoring, weather data and fungicide application records.

The scripts were mainly used to generate figures and summary tables for the Tiptree strawberry analysis chapter. They are not currently arranged as a single automated pipeline, so individual scripts should be run as needed after checking input paths and filenames.

## Repository structure

- `Supplementary_Data/`: supporting data files used by some of the analysis scripts.
- `*.R`: R scripts for data processing, plotting and statistical exploration.
- `*.py`: Python helper scripts for parsing alignment outputs.
- `*.ipynb`: exploratory Jupyter notebooks.

## Main scripts

### `Tiptree_all_data_analysis.R`

Main figure-generation script for the Tiptree analysis. It combines December sequencing data with the February to August PHI-base mapping results, focuses on the key pathogen genera `Botrytis`, `Podosphaera` and `Phytophthora`, and generates plots for:

- Normalised sequencing hits over time.
- Disease scores over time.
- Temperature and humidity summaries.
- Fungicide application timing.
- Combined multi-panel figures showing sequencing, disease, weather and fungicide data together.

### `PHIbase_analysis.R`

Analyses PHI-base mapping results across the full dataset. It normalises read counts, identifies the most abundant species and creates:

- Top 20 species bar charts.
- Bar charts with error bars.
- Stacked bar charts by month and location.
- A summary table of species counts under different filtering thresholds.

### `Heatmaps.R`

Creates heatmaps from PHI-base alignment data. The script filters taxa by read count and detection frequency, normalises hits per 100,000 reads, and produces heatmaps grouped by:

- Individual sample.
- Month.
- Location.
- Combined month and location plots.

### `sequencing_disease_merge.R`

Exploratory script for combining sequencing and disease score data. It merges normalised sequencing hits with disease scores by location, date and genus, then explores relationships between sequencing detections and disease scoring using:

- Summary statistics.
- Histograms.
- Scatter plots and trend lines.
- Monthly aggregation.
- Pearson correlation tests.
- Exploratory cross-correlation analysis.

### `fungicide_data.R`

Processes fungicide spray records and FRAC resistance-risk information. It summarises spray applications by location, month and target disease, then creates plots showing:

- Number of fungicide applications per month.
- Fungicide target categories.
- FRAC resistance-risk categories by location.

### `DNA yield Boxplot.R`

Creates boxplots of DNA yield by collection date. It compares first extraction yield and post-WGA yield, including a combined plot showing both yield measurements on the same graph.

### `Tiptree_DNA_yield_reads_barcode.R`

Explores the relationship between DNA yield and sequencing output. It creates scatter plots and linear model trend lines comparing read count with:

- Extraction DNA yield.
- Logged extraction DNA yield.
- T7/WGA yield.
- Logged T7/WGA yield.

### `old_Tiptree_Q20_analysis_Feb_script.R`

Older exploratory analysis script for the February to August Q20 PHI-base mapping data. It prepares genus-level summaries for `Botrytis`, `Podosphaera` and `Phytophthora`, writes out a processed CSV, and includes early versions of plots for sequencing hits, disease scores, fungicide applications, weather data and heatmaps.

### `error_bar_sequence_Richards_script.R`

Older reference-style plotting script for AirSeq pathogen hit data and weather data. It generates line plots with error bars for selected pathogen species and combines these plots into a multi-panel PDF.

## Python helper scripts

### `dict_paf_parse.py`

Parses a PAF alignment file and counts reads assigned to each taxa ID. It stores alignments by query name, filters on mapping quality, counts unique mappings directly, resolves some multi-mapping reads using mapping quality, and reports reads that were ignored because they mapped equally to multiple taxa.

### `SAM_counter.py`

Counts alignments from a SAM file by reference group. It uses a FASTQ file to calculate read lengths and a reference mapping file to connect sequence names to reference names, then outputs the number of aligned reads and total aligned base pairs per reference.

## Jupyter notebooks

### `data_statistics.ipynb`

Exploratory notebook for looking at sample metadata. It reads the Tiptree metadata table, separates air and mildew samples, cleans DNA yield fields, converts collection dates and generates descriptive statistics for DNA yield and sequencing metrics.

### `mildew_pipeline_anaysis.ipynb`

Exploratory notebook for checking mildew/emergent pathogen pipeline outputs. It examines genome coverage results, focuses on mildew-related taxa such as `Podosphaera aphanis`, and explores coverage summaries from combined mappings and `samtools depth` output.

## Requirements

The R scripts use packages including:

- `ggplot2`
- `dplyr`
- `tidyr`
- `lubridate`
- `scales`
- `patchwork`
- `RColorBrewer`
- `forcats`
- `stringr`
- `broom`
- `zoo`
- `astsa`
- `grid`
- `gridExtra`

The Python notebooks use packages including:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
