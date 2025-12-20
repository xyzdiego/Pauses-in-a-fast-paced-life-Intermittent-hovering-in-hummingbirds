# Pauses in a fast paced life intermittent hovering in hummingbirds

![Behavioral Ecophysics Banner](https://ecophysics.org/wp-content/uploads/2024/09/Behavioral-Ecophysics-text-logo-3D.png)

This repository contains the source code and data required to reproduce comparative phylogenetic analyses on behavioral (flight pauses) and morphological (wing coloration, morphometrics) traits in a hummingbird species assemblage (*Trochilidae*).

## Project Description

The goal of this project is to evaluate the evolutionary structure and adaptive correlations between flight behavior and morphology. The main script (`analysis.R`) executes a comprehensive workflow ranging from data cleaning to generating publication-ready figures.

Key analyses include:
* **Data & Phylogenetic Curation:** Taxonomic synchronization and manual addition of missing species to the phylogenetic tree using `bind.tip()`.
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

This analysis is executed in **R**. Key dependencies include standard CRAN packages and `ggtree` from Bioconductor. Run the following code to set up the environment:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("ggtree", quietly = TRUE)) BiocManager::install("ggtree")

if(!require(pacman)) install.packages("pacman")
pacman::p_load(phytools, ggtree, caper, tidyverse, RColorBrewer, 
               readxl, janitor, nlme, ape, phylolm, cowplot, 
               caret, kableExtra, ggplotify)
```

## Visualization & Outputs

The code automatically generates high-quality composite figures:
 
 * **Tanglegrams:** Mirrored visualization comparing the evolutionary history of flight pauses vs. underwing coloration.
 * **Composite Plots:** Integration of phylogenies with violin plots (data distribution) and forest plots (PGLS coefficients with 95% confidence intervals).
 * **Summary Tables:** Stylized tables containing comparative metrics (q rates, D statistic, F-statistics, AIC).

## Contribution & Usage

If you use this code or data for your research, please cite this repository and the associated article.

License: MIT
