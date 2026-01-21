# ==============================================================================
# AUXILIARY FUNCTIONS SCRIPT
# Project: Pauses in a fast paced life: Intermittent hovering in hummingbirds
# Description: Contains custom functions for phylogenetic models,
#              visualization (forest plots), and comparative metrics calculation.
#
# Required dependencies (must be loaded in the main script):
#   - tidyverse (dplyr, ggplot2)
#   - caper, phytools, phylolm
#
# ==============================================================================

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
    # Requires phylolm loaded
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
    # Requires ggplot2 loaded
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
    # WARNING: This function assumes the dataframe has columns named exactly:
    # "wing_length", "body_mass", "winglength_vs_mass"
    
    morpho_vars <- c("wing_length", "body_mass", "winglength_vs_mass")
    
    # Explicit use of dplyr:: to avoid conflicts with other packages like MASS or raster
    data_base <- df_full %>%
        dplyr::select(species, dplyr::all_of(col_name), dplyr::all_of(morpho_vars)) %>%
        dplyr::rename(trait_cat = dplyr::all_of(col_name)) %>%
        dplyr::filter(trait_cat %in% c("YES", "NO")) %>%
        tidyr::drop_na() %>%
        droplevels()
    
    data_base <- data.frame(data_base)
    sp_temp <- intersect(data_base$species, phy_tree$tip.label)
    tree_temp <- ape::keep.tip(phy_tree, sp_temp)
    data_base <- data_base[match(sp_temp, data_base$species), ]
    
    states_vec <- setNames(as.character(data_base$trait_cat), 
                           data_base$species)
    fit_mk <- phytools::fitMk(tree_temp, states_vec, model = "ER")
    
    data_base$trait_bin <- ifelse(data_base$trait_cat == "YES", 1, 0)
    comp_data_d <- caper::comparative.data(tree_temp, data_base, names.col = "species")
    d_stat <- caper::phylo.d(comp_data_d, binvar = trait_bin, permut = 200)
    
    results <- c(
        "N (Species)"      = nrow(data_base),
        "ACE Rate (q)"     = sprintf("%.4f", fit_mk$rates),
        "Phylo Signal (D)" = sprintf("%.3f", d_stat$DEstimate),
        "D P-value"        = sprintf("%.4f", d_stat$Pval1)
    )
    
    comp_data_pgls <- caper::comparative.data(tree_temp, data_base, 
                                              names.col = "species", vcv = TRUE)
    
    for (var in morpho_vars) {
        form_pgls <- as.formula(paste(var, "~ trait_cat"))
        mod_pgls <- caper::pgls(form_pgls, data = comp_data_pgls, lambda = "ML")
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