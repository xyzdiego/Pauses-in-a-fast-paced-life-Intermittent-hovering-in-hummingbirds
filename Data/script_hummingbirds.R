# SET WORKING DIRECTORY ---------------------------------------------------

WORK_DIR <- "HERE" # Please, insert the address of the directory that contains the data
setwd(WORK_DIR)

# LIBRARIES AND CONFIGURATION --------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(
    phytools, ggtree, caper, tidyverse, RColorBrewer, readxl, janitor, nlme, 
    ape, phylolm, cowplot, caret, kableExtra, ggplotify
)

# AUXILIARY FUNCTIONS ----------------------------------------------------------

# Function to extract coefficients from phyloglm models
extract_coefficients <- function(model, category) {
    if (is.null(model)) return(NULL)
    
    suma <- try(summary(model), silent = TRUE)
    if (inherits(suma, "try-error") || !"coefficients" %in% names(suma)) 
        return(NULL)
    
    coefs <- suma$coefficients
    if (is.null(coefs) || nrow(coefs) == 0) return(NULL)
    
    result <- data.frame(
        Category = rep(category, nrow(coefs)),
        Variable = rownames(coefs),
        stringsAsFactors = FALSE
    )
    
    n_cols <- ncol(coefs)
    if (n_cols >= 1) result$Coefficient <- round(coefs[, 1], 4)
    if (n_cols >= 2) result$Std_Error   <- round(coefs[, 2], 4)
    if (n_cols >= 3) result$Statistic   <- round(coefs[, 3], 4)
    if (n_cols >= 4) result$P_Value     <- round(coefs[, 4], 4)
    
    return(result)
}

# Function to generate data for Forest Plots (Confidence Intervals)
get_forest_data <- function(dep_var, pred_var, data, tree, colors) {
    levels_vec <- c("NO", "YES")
    res_df <- data.frame()
    
    for (cat in levels_vec) {
        data$y_bin <- ifelse(data[[dep_var]] == cat, 1, 0)
        form <- as.formula(paste("y_bin ~", pred_var))
        mod <- phyloglm(form, data = data, phy = tree, method = "logistic_MPLE")
        
        cfs <- summary(mod)$coefficients
        row_idx <- which(rownames(cfs) != "(Intercept)")[1]
        
        est <- cfs[row_idx, "Estimate"]
        se <- if ("StdErr" %in% colnames(cfs)) {
            cfs[row_idx, "StdErr"] 
        } else {
            cfs[row_idx, 2]
        }
        
        tmp <- data.frame(
            Category = cat, Variable = "Wing length",
            Estimate = est, CI_Low = est - 1.96 * se, CI_High = est + 1.96 * se
        )
        res_df <- rbind(res_df, tmp)
    }
    return(res_df)
}

# Function to create Forest Plots
create_forest_plot <- function(df, cols, title_txt) {
    ggplot(df, aes(x = Estimate, y = Variable, color = Category)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_errorbar(aes(xmin = CI_Low, xmax = CI_High), 
                      width = 0.2, size = 1, 
                      position = position_dodge(width = 0.5)) +
        geom_point(size = 5, position = position_dodge(width = 0.5)) +
        scale_color_manual(values = cols) +
        labs(title = title_txt, x = "Confidence interval ±95%", y = " ") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              axis.text.y = element_text(size = 12, face = "bold", 
                                         angle = 90, hjust = 0.5),
              legend.position = "right",
              panel.grid.major.y = element_blank())
}

# Function to calculate comparative metrics for the summary table
calculate_full_metrics <- function(col_name, df_full, phy_tree) {
    morpho_vars <- c("wing_length", "body_mass", "winglength_vs_mass")
    
    data_base <- df_full %>%
        dplyr::select(species, all_of(col_name), all_of(morpho_vars)) %>%
        rename(trait_cat = all_of(col_name)) %>%
        filter(trait_cat %in% c("YES", "NO")) %>%
        drop_na() %>%
        droplevels()
    
    data_base <- data.frame(data_base)
    sp_temp <- intersect(data_base$species, phy_tree$tip.label)
    tree_temp <- keep.tip(phy_tree, sp_temp)
    data_base <- data_base[match(sp_temp, data_base$species), ]
    
    states_vec <- setNames(as.character(data_base$trait_cat), 
                           data_base$species)
    fit_mk <- fitMk(tree_temp, states_vec, model = "ER")
    
    data_base$trait_bin <- ifelse(data_base$trait_cat == "YES", 1, 0)
    comp_data_d <- comparative.data(tree_temp, data_base, names.col = "species")
    d_stat <- phylo.d(comp_data_d, binvar = trait_bin, permut = 200)
    
    results <- c(
        "N (Species)"      = nrow(data_base),
        "ACE Rate (q)"     = sprintf("%.4f", fit_mk$rates),
        "Phylo Signal (D)" = sprintf("%.3f", d_stat$DEstimate),
        "D P-value"        = sprintf("%.4f", d_stat$Pval1)
    )
    
    comp_data_pgls <- comparative.data(tree_temp, data_base, 
                                       names.col = "species", vcv = TRUE)
    
    for (var in morpho_vars) {
        form_pgls <- as.formula(paste(var, "~ trait_cat"))
        mod_pgls <- pgls(form_pgls, data = comp_data_pgls, lambda = "ML")
        sum_pgls <- summary(mod_pgls)
        
        f_stat <- sum_pgls$fstatistic[1]
        p_val <- pf(f_stat, sum_pgls$fstatistic[2], sum_pgls$fstatistic[3], 
                    lower.tail = FALSE)
        
        new_res <- c(
            "PGLS Lambda" = sprintf("%.3f", sum_pgls$param["lambda"]),
            "PGLS F-stat" = sprintf("%.2f", f_stat),
            "PGLS P-value" = sprintf("%.5f", p_val),
            "PGLS AIC"    = sprintf("%.1f", AIC(mod_pgls))
        )
        results <- c(results, new_res)
    }
    return(results)
}

# DATA IMPORT AND PROCESSING ---------------------------------------------------

data_raw <- read_excel("Database_Simulation.xlsx", sheet = "Pausas_Tipo")
data_raw <- data_raw[data_raw$Species != "Sternoclyta cyanopectus", ]
tree_raw <- read.tree("arbol.nwk") 

data_raw[, c(FALSE, sapply(data_raw[,-1], is.character))] <- 
    lapply(data_raw[, c(FALSE, sapply(data_raw[,-1], is.character))], factor)
data_clean <- data_raw %>% clean_names()

data_clean$species <- gsub(" ", "_", data_clean$species)
tree_raw$tip.label <- gsub(" ", "_", tree_raw$tip.label)

tree_mod <- bind.tip(tree_raw, "Aglaeactis_aliciae", edge.length = 0.2, 
                     where = 176)
tree_mod <- bind.tip(tree_mod, "Lamprolaima_rhami", edge.length = 3.2, 
                     where = 107)
tree_mod <- bind.tip(tree_mod, "Sternoclyta_cyanopectus", edge.length = 6.0, 
                     where = 108)
tree_mod <- bind.tip(tree_mod, "Oxypogon_lindenii", edge.length = 1.6, 
                     where = 230)
tree_mod <- bind.tip(tree_mod, "Urochroa_leucura", edge.length = 0.8, 
                     where = 154)

name_replacements <- c(
    "Chlorostilbon_notatus" = "Chlorestes_notata", 
    "Damophila_julie" = "Chlorestes_julie",
    "Goethalsia_bella" = "Goldmania_bella", 
    "Lophornis_pavoninus" = "Lophornis_adorabilis",
    "Cyanophaia_bicolor" = "Riccordia_bicolor", 
    "Amazilia_beryllina" = "Saucerottia_beryllina",
    "Amazilia_castaneiventris" = "Saucerottia_castaneiventris", 
    "Amazilia_cyanura" = "Saucerottia_cyanura",
    "Heliodoxa_cyanopectus" = "Sternoclyta_cyanopectus", 
    "Campylopterus_villaviscencio" = "Campylopterus_villaviscensio"
)

tree_mod$tip.label <- ifelse(tree_mod$tip.label %in% names(name_replacements), 
                             name_replacements[tree_mod$tip.label], 
                             tree_mod$tip.label)

if (any(tree_mod$edge.length <= 0)) {
    tree_mod$edge.length[tree_mod$edge.length <= 0] <- 1e-6
}
if (!is.ultrametric(tree_mod)) tree_mod <- chronopl(tree_mod, lambda = 1)
tree_mod$tip.label <- gsub(" ", "_", tree_mod$tip.label)

common_spp <- intersect(data_clean$species, tree_mod$tip.label)
tree_bin <- keep.tip(tree_mod, common_spp)
data_bin <- data_clean %>% 
    filter(species %in% common_spp) %>% 
    arrange(match(species, tree_bin$tip.label))

data_bin <- as.data.frame(data_bin)
rownames(data_bin) <- data_bin$species

data_bin <- data_bin %>%
    filter(hovering_pauses %in% c("YES", "NO")) %>%
    filter(color %in% c("YES", "NO")) %>%
    filter(!is.na(wing_length)) %>% 
    mutate(winglength_vs_mass = winglength_rel_mass2)

final_spp <- intersect(data_bin$species, tree_bin$tip.label)
tree_bin <- keep.tip(tree_bin, final_spp)
data_bin <- data_bin[match(tree_bin$tip.label, data_bin$species), ]

if (any(duplicated(tree_bin$tip.label))) {
    dup_idx <- which(duplicated(tree_bin$tip.label))
    tree_bin <- drop.tip(tree_bin, dup_idx)
    data_bin <- data_bin[match(tree_bin$tip.label, data_bin$species), ]
}

data_bin$hovering_pauses_01 <- ifelse(data_bin$hovering_pauses == "YES", 1, 0)

levels_bin <- c("NO", "YES")
cols_bin <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(cols_bin) <- levels_bin

# PHYLOGENETIC SIGNAL ANALYSIS (D) ---------------------------------------------

signals_bin <- list()

for (cat in levels_bin) {
    var_bin <- as.numeric(data_bin$hovering_pauses == cat)
    names(var_bin) <- data_bin$species
    
    data_temp <- data.frame(
        species = data_bin$species, 
        state_bin = var_bin, 
        stringsAsFactors = FALSE
    )
    
    comp_data <- comparative.data(tree_bin, data_temp, names.col = "species")
    d_res <- phylo.d(comp_data, binvar = state_bin, permut = 100)
    signals_bin[[cat]] <- d_res
}

# PHYLOGLM MODELING (HOVERING) -------------------------------------------------

model_hover_yes <- phyloglm(hovering_pauses_01 ~ wing_length, 
                            data = data_bin, phy = tree_bin, 
                            method = "logistic_MPLE")

# Coefficient Extraction
model_list <- list("Hovering (YES)" = model_hover_yes)
all_coefs <- data.frame()
sig_effects <- data.frame()

for (cat in names(model_list)) {
    cat_coefs <- extract_coefficients(model_list[[cat]], cat)
    if (!is.null(cat_coefs) && nrow(cat_coefs) > 0) {
        all_coefs <- rbind(all_coefs, cat_coefs)
        if ("P_Value" %in% colnames(cat_coefs)) {
            sig <- cat_coefs[cat_coefs$P_Value < 0.05, ]
            if (nrow(sig) > 0) sig_effects <- rbind(sig_effects, sig)
        }
    }
}

# PROBABILITY HEATMAP ----------------------------------------------------------

prob_yes <- model_hover_yes$fitted.values
prob_matrix <- cbind(NO = 1 - prob_yes, YES = prob_yes)
rownames(prob_matrix) <- data_bin$species

if (all(!is.na(prob_matrix))) {
    phylo.heatmap(tree_bin, prob_matrix, standardize = FALSE, 
                  legend = TRUE, fsize = 1.2, colors = c("#E41A1C", "#377EB8"))
    title("Predicted Probabilities (P(NO), P(YES))")
}

# EVOLUTIONARY CORRELATION (FITPAGEL) ------------------------------------------

# A. Correlation vs Microhabitat
data_pagel_hab <- data_bin %>%
    mutate(habitat_bin = case_when(
        microhabitat %in% c("Canopy", "Understory") ~ "Closed",
        microhabitat %in% c("Open", "Mixed")        ~ "Open",
        TRUE ~ NA_character_ 
    )) %>%
    filter(!is.na(habitat_bin), !is.na(hovering_pauses)) %>%
    droplevels()

sp_hab <- intersect(data_pagel_hab$species, tree_bin$tip.label)
tree_hab <- keep.tip(tree_bin, sp_hab)
data_pagel_hab <- data_pagel_hab[match(sp_hab, data_pagel_hab$species), ]

r1 <- setNames(factor(data_pagel_hab$hovering_pauses), data_pagel_hab$species)
r2 <- setNames(factor(data_pagel_hab$habitat_bin), data_pagel_hab$species)

if (length(levels(r1)) == 2 && length(levels(r2)) == 2){
    fit_dep_hab <- fitPagel(tree_hab, r1, r2, model = "ER")
}

# B. Correlation vs Foraging
data_pagel_for <- data_bin %>%
    mutate(forage_bin = case_when(
        foraging_type == "Territorial" ~ "Territorial",
        foraging_type == "Opportunist" ~ "No_Territorial",
        TRUE ~ NA_character_ 
    )) %>%
    filter(!is.na(forage_bin), !is.na(hovering_pauses)) %>%
    droplevels()

sp_for <- intersect(data_pagel_for$species, tree_bin$tip.label)
tree_for <- keep.tip(tree_bin, sp_for)
data_pagel_for <- data_pagel_for[match(sp_for, data_pagel_for$species), ]

r3 <- setNames(factor(data_pagel_for$forage_bin), data_pagel_for$species)

if (length(levels(r1)) == 2 && length(levels(r3)) == 2){
    fit_dep_for <- fitPagel(tree_for, r1, r3, model = "ER")
}

# TANGLEGRAM PLOT (MIRRORED TREE) ----------------------------------------------

# Colors and Clades
unique_clades <- unique(data_bin$clade)
cols_clades <- RColorBrewer::brewer.pal(max(length(unique_clades), 3), 
                                        "Paired")[1:length(unique_clades)]
names(cols_clades) <- unique_clades

edge_cols <- rep("gray10", nrow(tree_bin$edge))
edge_widths <- rep(1, nrow(tree_bin$edge))

for (c in unique_clades) {
    sp_clade <- intersect(data_bin$species[data_bin$clade == c], 
                          tree_bin$tip.label)
    if (length(sp_clade) > 0) {
        idx <- which.edge(tree_bin, sp_clade)
        edge_cols[idx] <- cols_clades[c]
        edge_widths[idx] <- 4
    }
}

# Wing Color Mapping (Corrected and aligned)
wing_levels <- unique(trimws(tolower(as.character(data_bin$wing_color))))
cols_wing_map <- setNames(wing_levels, wing_levels)

if ("indigo" %in% names(cols_wing_map)) cols_wing_map["indigo"] <- "#4B0082"
if ("sapphire" %in% names(cols_wing_map)) cols_wing_map["sapphire"] <- "#2138AB"
if ("white" %in% names(cols_wing_map)) cols_wing_map["white"] <- "gray60"
if ("gray" %in% names(cols_wing_map)) cols_wing_map["gray"] <- "gray35"
if ("orange" %in% names(cols_wing_map)) cols_wing_map["orange"] <- "darkorange"

# Explicit vector extraction to avoid tibble indexing errors
ordered_idx <- match(tree_bin$tip.label, data_bin$species)
ordered_colors <- as.character(data_bin$wing_color[ordered_idx])
cols_wing_sp <- cols_wing_map[trimws(tolower(ordered_colors))]

# Plot Setup
max_depth <- max(node.depth.edgelength(tree_bin))
xlim_l <- c(-max_depth * 0.35, max_depth * 2.6) 
xlim_r <- c(-max_depth * 0.35, max_depth * 2.6)

clean_labels <- gsub("_", " ", tree_bin$tip.label)
clean_labels <- gsub("Coeligena bonapartei", "Coeligena bonapartei/consita", 
                     clean_labels)

layout(matrix(c(1,2), 1, 2))
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 2, 0)) 

# Left: Hovering Pauses
states_hover <- setNames(as.character(data_bin$hovering_pauses), 
                         data_bin$species)
ace_hover <- ace(states_hover, tree_bin, type = "discrete", model = "ER")
cols_hover_fix <- c("NO" = "#1B9E77", "YES" = "#D95F02") 

plot(tree_bin, cex = 0.95, no.margin = TRUE, label.offset = 0.02, 
     direction = "rightwards", show.tip.label = FALSE, 
     x.lim = xlim_l, edge.color = edge_cols, edge.width = edge_widths)

mtext("Hovering Pauses", side = 3, line = -2, font = 2, at = max_depth * 0.5)
nodelabels(pie = ace_hover$lik.anc, piecol = cols_hover_fix, cex = 0.45)
tiplabels(pch = 21, bg = cols_hover_fix[states_hover], cex = 3)
legend("bottomleft", legend = names(cols_hover_fix), fill = cols_hover_fix,
       title = "Hovering", cex = 1, bty = "n")

# Center: Labels
pos_center <- max_depth * (2.6 - 1) * 0.5 + (max_depth * 0.1)
tiplabels(text = clean_labels, frame = "none", col = cols_wing_sp, 
          font = 2, cex = 0.9, adj = 0.5, offset = pos_center)

# Right: Underwing Coloration
if (!exists("ace_color")) {
    states_col <- setNames(as.character(data_bin$color), data_bin$species)
    ace_color <- ace(states_col, tree_bin, type = "discrete", model = "ER")
}

plot(tree_bin, cex = 0.95, no.margin = TRUE, label.offset = 0.02, 
     direction = "leftwards", show.tip.label = FALSE, 
     x.lim = xlim_r, edge.color = edge_cols, edge.width = edge_widths)

mtext("Underwing coloration", side = 3, line = -2, font = 2, at = max_depth * 0.5)
cols_right <- c("NO" = "gray40", "YES" = "#4B0082")

nodelabels(pie = ace_color$lik.anc, piecol = cols_right, cex = 0.45)
tiplabels(pch = 21, bg = cols_right[data_bin$color], cex = 3)
legend("bottomright", legend = c("Gray", "Colorful"), 
       fill = c("gray40", "#4B0082"), title = "Underwing coloration", 
       cex = 1, bty = "n")

legend("topright", legend = names(cols_clades), col = cols_clades, 
       lwd = 5, title = "Clades", cex = 0.8, bty = "n", ncol = 1)

par(mfrow = c(1,1))

# COMPOSITE PLOT (VIOLINS + CONFIDENCE INTERVALS) ------------------------------

# Colors
cols_hover_plt <- RColorBrewer::brewer.pal(3, "Dark2")[1:2]
names(cols_hover_plt) <- c("NO", "YES")
cols_color_plt <- c("NO" = "gray40", "YES" = "#4B0082")

# Data preparation
plot_data <- data_bin %>%
    filter(hovering_pauses %in% c("YES", "NO"),
           color %in% c("YES", "NO"),
           !is.na(wing_length)) %>%
    droplevels()

# Violin Plots
p_violin_hover <- ggplot(plot_data, aes(x = hovering_pauses, y = wing_length, 
                                        fill = hovering_pauses)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.4) +
    scale_fill_manual(values = cols_hover_plt) +
    labs(title = "Wing Length vs Hovering Pauses",
         y = "Wing Length (mm)", x = "Hovering Pauses") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

p_violin_color <- ggplot(plot_data, aes(x = color, y = wing_length, 
                                        fill = color)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.4) +
    scale_fill_manual(values = cols_color_plt) +
    scale_x_discrete(labels = c("NO" = "Gray", "YES" = "Colorful")) +
    labs(title = "Wing Length vs Underwing Coloration",
         y = "Wing Length (mm)", x = "Underwing Coloration") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

# Forest Plots generation
df_conf_hover <- get_forest_data("hovering_pauses", "wing_length", data_bin, 
                                 tree_bin, cols_hover_plt)
df_conf_color <- get_forest_data("color", "wing_length", data_bin, 
                                 tree_bin, cols_color_plt)

p_conf_hover <- create_forest_plot(df_conf_hover, cols_hover_plt, 
                                   "Hovering Pauses")
p_conf_color <- create_forest_plot(df_conf_color, cols_color_plt, 
                                   "Underwing Coloration")

final_plot <- plot_grid(p_conf_hover, p_violin_hover, p_conf_color, 
                        p_violin_color, ncol = 2, align = 'h', 
                        labels = c("A", "B", "C", "D"))
print(final_plot)

# COMPREHENSIVE SUMMARY TABLE --------------------------------------------------

res_hover <- calculate_full_metrics("hovering_pauses", data_bin, tree_bin)
res_trans <- calculate_full_metrics("transition_pauses", data_bin, tree_bin)
res_color <- calculate_full_metrics("color", data_bin, tree_bin)

summary_table <- data.frame(
    Metric = names(res_hover),
    Hovering = res_hover,
    Transitions = res_trans,
    Coloration = res_color,
    row.names = NULL,
    stringsAsFactors = FALSE
)

# Highlight logic: All 3 models must be significant (< 0.05)
idx_pvals <- c(7, 11, 15)
red_rows <- c()

for (i in idx_pvals) {
    v_h <- as.numeric(summary_table[i, "Hovering"])
    v_t <- as.numeric(summary_table[i, "Transitions"])
    v_c <- as.numeric(summary_table[i, "Coloration"])
    if (v_h < 0.05 & v_t < 0.05 & v_c < 0.05) {
        red_rows <- c(red_rows, i)
    }
}

kbl_table <- summary_table %>%
    kbl(caption = "Comparative Summary of Evolutionary Statistics and PGLS Models", 
        align = "c") %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    row_spec(0, bold = T, color = "white", background = "#2c3e50") %>%
    column_spec(1, bold = T, border_right = T, width = "4cm", color = "#555") %>%
    pack_rows("1. Evolutionary Signal (Categorical Trait)", 1, 4, 
              label_row_css = "background-color:#ecf0f1;color:#333;font-weight:bold;") %>%
    pack_rows("2. Model: Wing Length ~ Trait", 5, 8, 
              label_row_css = "background-color:#d1f2eb;color:#0e6655;font-weight:bold;") %>%
    pack_rows("3. Model: Body Mass ~ Trait", 9, 12, 
              label_row_css = "background-color:#d6eaf8;color:#154360;font-weight:bold;") %>%
    pack_rows("4. Model: Relation Wing Length vs Body Mass ~ Trait", 13, 16, 
              label_row_css = "background-color:#fce8d6;color:#a04000;font-weight:bold;") %>%
    row_spec(idx_pvals, bold = T)

if (length(red_rows) > 0) {
    kbl_table <- kbl_table %>% row_spec(red_rows, color = "red", bold = T)
}

print(kbl_table)

# Transtition pauses ------------------------------------------------------

cols_trans <- c("NO" = "#1B9E77", "YES" = "#D95F02") 

data_trans_plot <- data_bin %>%
    filter(transition_pauses %in% c("YES", "NO"),
           !is.na(wing_length)) %>%
    droplevels()

states_trans <- setNames(as.character(data_trans_plot$transition_pauses), 
                         data_trans_plot$species)
tree_trans <- keep.tip(tree_bin, names(states_trans))
ace_trans <- ace(states_trans, tree_trans, type = "discrete", model = "ER")

p_tree_trans <- as.grob(function() {
    par(mar = c(0, 0, 2, 0)) 
    plot(tree_trans, 
         cex = 0.8,                   
         no.margin = TRUE, 
         label.offset = 0.05, 
         direction = "rightwards", 
         show.tip.label = TRUE,       
         edge.color = edge_cols,      
         edge.width = 3)              
    mtext("Transition Pauses Phylogeny", side = 3, line = -1, font = 2, cex = 1)
    nodelabels(pie = ace_trans$lik.anc, piecol = cols_trans, cex = 0.4)
    tiplabels(pch = 21, bg = cols_trans[states_trans], cex = 1.2, offset = 0.02)
    legend("bottomleft", legend = names(cols_trans), fill = cols_trans,
           title = "Transitions", cex = 0.9, bty = "n")
})

p_violin_trans <- ggplot(data_trans_plot, aes(x = transition_pauses, 
                                              y = wing_length, 
                                              fill = transition_pauses)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.4) +
    scale_fill_manual(values = cols_trans) +
    labs(title = "Wing Length Distribution",
         y = "Wing Length (mm)", x = "Transition Pauses") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

df_conf_trans <- get_forest_data("transition_pauses", "wing_length", 
                                 data_trans_plot, tree_trans, cols_trans)
p_forest_trans <- create_forest_plot(df_conf_trans, cols_trans, 
                                     "Evolutionary Effect (Model)")

col_right <- plot_grid(p_violin_trans, p_forest_trans, 
                         ncol = 1, 
                         labels = c("B", "C"))

final_trans_plot <- plot_grid(p_tree_trans, col_right, 
                              ncol = 2, 
                              labels = c("A", ""), 
                              rel_widths = c(1, 1)) 

print(final_trans_plot)
