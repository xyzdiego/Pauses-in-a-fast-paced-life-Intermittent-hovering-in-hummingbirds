# Data Dictionary: Trochilidae Morphology & Behavior

**Dataset File:** `Trochilidae_MorphoBehavior_Data.xlsx` (Sheet: "Pausas_Tipo")
**Description:** Comparative dataset containing morphological measurements, flight mechanics, and behavioral classifications for hummingbird species.

## Variable Descriptions

| Variable Name (Raw) | Data Type | Units | Description | Possible Values / Notes |
| :--- | :--- | :--- | :--- | :--- |
| `Species` | String | - | Scientific name of the species | Format: *Genus species* |
| `Hovering pauses` | Categorical | - | Presence of intermittent hovering behavior | **YES**, **NO**, **U** (Unknown/Uncertain) |
| `Transition pauses`| Categorical | - | Presence of pauses during flight transitions | **YES**, **NO**, **U** (Unknown/Uncertain) |
| `Body mass` | Numeric | g | Total body mass | |
| `Wing area` | Numeric | mm² | Total area of both wings | Measurement with large variance for the type of units |
| `Wing loading` | Numeric | g/cm² | Ratio of body mass to wing area | |
| `wingLength` | Numeric | mm | Length of the wing chord | |
| `TailLength` | Numeric | mm | Length of the tail | |
| `totalLength` | Numeric | mm | Total length of the bird (bill to tail) | |
| `BodyLength` | Numeric | mm | Length of the body (excluding bill/tail) | |
| `Color` | Categorical | - | Presence of iridescent underwing coloration | **YES**: Colorful; **NO**: Drab/Gray |
| `Colorful` | Numeric | Score | Quantified score of coloration intensity | Scale 1-3 (Estimated) |
| `colorful %` | Numeric | % | Percentage of wing area covered by coloration | 0 - 100 |
| `Wing color` | Categorical | - | Dominant color description of the wing | e.g., "Gray", "White", "Orange" |
| `Tail color` | Categorical | - | Dominant color description of the tail | e.g., "Brown", "Black" |
| `MaxAlt` | Numeric | m | Maximum elevation of species distribution | Meters above sea level |
| `MinAlt` | Numeric | m | Minimum elevation of species distribution | Meters above sea level |
| `RangeAlt` | Numeric | m | Elevational range (Max - Min) | |
| `Clade` | Categorical | - | Taxonomic clade assignment | e.g., "Coquettes", "Brilliants" |
| `Foraging type` | Categorical | - | Foraging strategy classification | "Traplining", "Territorial", "Generalist" |
| `Microhabitat` | Categorical | - | Preferred forest stratum | "Understory", "Mixed", "Canopy" |
| `Muscle type` | Numeric | Type | Classification of muscle fiber composition | Encoded as integer (e.g., 2, 3) |
| `winglength_relMass2`| Numeric | Index | Size-corrected wing index (Model 2) | Residuals of Wing Length ~ Body Mass |
| `winglength_relMass3`| Numeric | Index | Size-corrected wing index (Model 3) | Alternative residual calculation |

## Data Processing Notes

1.  **Variable Renaming in R:** The analysis script applies the `clean_names()` function (from `janitor` package). This converts the raw column names listed above to `snake_case` for consistency in coding.
    * *Example:* `Hovering pauses` becomes `hovering_pauses`.
    * *Example:* `wingLength` becomes `wing_length`.
2.  **Missing Data:** Cells containing `NA` represent missing measurements.
3.  **Filtering:** Values marked as "U" (Unknown) in behavioral traits are filtered out prior to PGLS modeling.
