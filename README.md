# Pauses-in-a-fast-paced-life-Intermittent-hovering-in-hummingbirds

![Behavioral Ecophysics Banner](https://ecophysics.org/wp-content/uploads/2024/09/Behavioral-Ecophysics-text-logo-3D.png)

This repository contains the source code and data required to reproduce comparative phylogenetic analyses on behavioral (flight pauses) and morphological (wing coloration, morphometrics) traits in a hummingbird species assemblage (*Trochilidae*).

## Project Description

The goal of this project is to evaluate the evolutionary structure and adaptive correlations between flight behavior and morphology. The main script (`analysis.R`) executes a comprehensive workflow ranging from data cleaning to generating publication-ready figures.

Key analyses include:
* **Data & Phylogenetic Curation:** Taxonomic synchronization and manual addition of missing species to the phylogenetic tree using `bind.tip`.
* **Ancestral State Reconstruction (ACE):** Inference of the evolutionary history of binary traits (presence/absence of pauses and coloration) using Maximum Likelihood (ML) models.
* **Phylogenetic Signal:** Evaluation of trait conservation using Fritz & Purvis' *D* statistic and Pagel's *Lambda*.
* **PGLS Models:** Fitting of Phylogenetic Generalized Least Squares models to test the association between morphology (wing length, body mass) and behavior, statistically controlling for evolutionary history.

## Repository Structure

* `Data/`:
  - `script_hummingbirds.R`: Main R script containing the full workflow (cleaning, modeling, plotting).
  - `Database_Simulation.xlsx`: Dataset containing morphological and behavioral traits.
  - `arbol.nwk`: Time-calibrated phylogeny in Newick format.
* `outputs/`: Directory for saving generated figures.

## Requirements & Installation

This analysis is executed in **R**. It requires the installation of the following libraries, which are managed via `pacman` within the script:

```r
install.packages("pacman")
pacman::p_load(phytools, ggtree, caper, tidyverse, RColorBrewer, 
               readxl, janitor, nlme, ape, phylolm, cowplot, 
               caret, kableExtra, ggplotify)
