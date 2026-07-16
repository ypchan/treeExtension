#' Check monophyly of taxa at a given metadata rank on a rooted tree
#'
#' @description
#' Assess whether taxa at a specified metadata rank (for example: phylum, class,
#' order, family, genus, species) are monophyletic on a rooted phylogenetic tree.
#'
#' The function aligns `data` to `tree$tip.label` using `tip_col`, then evaluates
#' each taxon at the requested `rank_col`.
#'
#' Classification logic:
#' - `monophyletic`: all members form one pure clade
#' - `near_monophyletic_minor_outliers`: one dominant cluster plus only a few
#'   stray members elsewhere
#' - `polyphyletic`: the taxon is split into multiple substantial clusters
#' - `non_monophyletic_split`: non-monophyletic, but not strongly polyphyletic
#' - `singleton` / `too_few_tips`: too few members to assess robustly
#'
#' @param tree A rooted `phylo` object.
#' @param data A data.frame or tibble containing tip metadata.
#' @param rank_col Character scalar. Column name in `data` giving the taxonomic
#'   rank to evaluate.
#' @param tip_col Character scalar. Column in `data` containing tip labels
#'   matching `tree$tip.label`. Default: `"genome_label"`.
#' @param taxa Optional character vector of taxa to evaluate. Default `NULL`
#'   means all taxa present at `rank_col`.
#' @param min_tips Minimum number of tips required to assess a taxon.
#'   Default: `2L`.
#' @param near_mono_frac Minimum fraction of members in the largest pure cluster
#'   to classify a taxon as `near_monophyletic_minor_outliers`. Default: `0.80`.
#' @param polyphyly_second_frac Minimum fraction of members in the second-largest
#'   pure cluster to classify a taxon as `polyphyletic`. Default: `0.30`.
#' @param polyphyly_largest_max Maximum fraction of members in the largest pure
#'   cluster for a taxon to be classified as `polyphyletic`. Default: `0.60`.
#' @param drop_na Logical. If `TRUE`, taxa with missing values in `rank_col`
#'   are ignored. Default: `TRUE`.
#' @param check_rooted Logical. If `TRUE`, error when `tree` is not rooted.
#'   Default: `TRUE`.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'   Default: `FALSE`.
#'
#' @return
#' A tibble with one row per taxon and these columns:
#' - `rank`
#' - `taxon`
#' - `n_tips`
#' - `status`
#' - `largest_pure_cluster_n`
#' - `largest_pure_cluster_frac`
#' - `second_pure_cluster_n`
#' - `second_pure_cluster_frac`
#' - `n_pure_clusters`
#' - `n_outside_largest_cluster`
#' - `mrca_clade_size`
#' - `mrca_contaminant_tips`
#' - `mrca_purity`
#' - `interpretation`
#'
#' @examples
#' # res_order <- check_rank_monophyly(
#' #   tree = tre,
#' #   data = metadata,
#' #   rank_col = "order",
#' #   tip_col = "genome_label"
#' # )
#' #
#' # dplyr::count(res_order, status)
#' # dplyr::filter(res_order, status != "monophyletic")
#'
#' @export
check_rank_monophyly <- function(
        tree,
        data,
        rank_col,
        tip_col = "genome_label",
        taxa = NULL,
        min_tips = 2L,
        near_mono_frac = 0.80,
        polyphyly_second_frac = 0.30,
        polyphyly_largest_max = 0.60,
        drop_na = TRUE,
        check_rooted = TRUE,
        verbose = FALSE
) {
    # ---- 1) Validate input ---------------------------------------------------
    if (!inherits(tree, "phylo")) {
        stop("`tree` must be a phylo object.")
    }
    
    if (!is.data.frame(data)) {
        stop("`data` must be a data.frame or tibble.")
    }
    
    if (!is.character(rank_col) || length(rank_col) != 1L || is.na(rank_col)) {
        stop("`rank_col` must be a single non-missing character string.")
    }
    
    if (!is.character(tip_col) || length(tip_col) != 1L || is.na(tip_col)) {
        stop("`tip_col` must be a single non-missing character string.")
    }
    
    if (!rank_col %in% colnames(data)) {
        stop("Column `", rank_col, "` was not found in `data`.")
    }
    
    if (!tip_col %in% colnames(data)) {
        stop("Column `", tip_col, "` was not found in `data`.")
    }
    
    if (!is.numeric(min_tips) || length(min_tips) != 1L || is.na(min_tips) || min_tips < 1) {
        stop("`min_tips` must be a single numeric value >= 1.")
    }
    min_tips <- as.integer(min_tips)
    
    if (!is.numeric(near_mono_frac) || length(near_mono_frac) != 1L ||
        is.na(near_mono_frac) || near_mono_frac < 0 || near_mono_frac > 1) {
        stop("`near_mono_frac` must be a single number between 0 and 1.")
    }
    
    if (!is.numeric(polyphyly_second_frac) || length(polyphyly_second_frac) != 1L ||
        is.na(polyphyly_second_frac) || polyphyly_second_frac < 0 || polyphyly_second_frac > 1) {
        stop("`polyphyly_second_frac` must be a single number between 0 and 1.")
    }
    
    if (!is.numeric(polyphyly_largest_max) || length(polyphyly_largest_max) != 1L ||
        is.na(polyphyly_largest_max) || polyphyly_largest_max < 0 || polyphyly_largest_max > 1) {
        stop("`polyphyly_largest_max` must be a single number between 0 and 1.")
    }
    
    if (!is.logical(drop_na) || length(drop_na) != 1L || is.na(drop_na)) {
        stop("`drop_na` must be TRUE or FALSE.")
    }
    
    if (!is.logical(check_rooted) || length(check_rooted) != 1L || is.na(check_rooted)) {
        stop("`check_rooted` must be TRUE or FALSE.")
    }
    
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.")
    }
    
    if (check_rooted && !ape::is.rooted(tree)) {
        stop("The tree is not rooted. Please root the tree before running this function.")
    }
    
    tip_ids <- trimws(as.character(data[[tip_col]]))
    if (any(is.na(tip_ids) | tip_ids == "")) {
        stop("Column `", tip_col, "` contains missing or empty tip IDs.")
    }
    
    if (anyDuplicated(tip_ids) > 0) {
        dup_ids <- unique(tip_ids[duplicated(tip_ids)])
        stop(
            "Duplicated tip IDs found in `", tip_col, "`. Example duplicates: ",
            paste(head(dup_ids, 5), collapse = ", ")
        )
    }
    
    # ---- 2) Align metadata to tree tip order --------------------------------
    tip_match <- match(tree$tip.label, tip_ids)
    missing_meta <- is.na(tip_match)
    
    if (any(missing_meta)) {
        stop(
            "Some tree tips are missing in `data`. Example missing tips: ",
            paste(head(tree$tip.label[missing_meta], 5), collapse = ", ")
        )
    }
    
    data_aligned <- data[tip_match, , drop = FALSE]
    rank_vec <- trimws(as.character(data_aligned[[rank_col]]))
    
    # ---- 3) Apply optional NA filtering --------------------------------------
    if (drop_na) {
        keep <- !is.na(rank_vec) & rank_vec != ""
        usable_tip_labels <- tree$tip.label[keep]
    } else {
        keep <- rep(TRUE, length(tree$tip.label))
        usable_tip_labels <- tree$tip.label
    }
    
    if (length(usable_tip_labels) < 2L) {
        stop("Fewer than 2 usable tips remain after filtering.")
    }
    
    # ---- 4) Reorder tree postorder for bottom-up counting --------------------
    tree2 <- ape::reorder.phylo(tree, order = "postorder")
    ntip <- length(tree2$tip.label)
    total_nodes <- ntip + tree2$Nnode
    edge <- tree2$edge
    
    # Re-align data to reordered tree
    tip_match2 <- match(tree2$tip.label, trimws(as.character(data_aligned[[tip_col]])))
    data_reordered <- data_aligned[tip_match2, , drop = FALSE]
    rank_vec2 <- trimws(as.character(data_reordered[[rank_col]]))
    
    if (drop_na) {
        usable_tip_idx <- which(!is.na(rank_vec2) & rank_vec2 != "")
    } else {
        usable_tip_idx <- seq_len(ntip)
    }
    
    # ---- 5) Build parent vector and descendant tip counts --------------------
    parent <- rep(NA_integer_, total_nodes)
    parent[edge[, 2]] <- edge[, 1]
    
    desc_total <- integer(total_nodes)
    desc_total[seq_len(ntip)] <- 1L
    
    for (i in seq_len(nrow(edge))) {
        p <- edge[i, 1]
        ch <- edge[i, 2]
        desc_total[p] <- desc_total[p] + desc_total[ch]
    }
    
    # ---- 6) Build mapping: taxon -> tip indices ------------------------------
    valid_rank <- rank_vec2
    
    if (drop_na) {
        valid_rank[!(seq_len(ntip) %in% usable_tip_idx)] <- NA
    }
    
    valid_taxon <- !is.na(valid_rank) & valid_rank != ""
    taxon2tips <- split(seq_len(ntip)[valid_taxon], valid_rank[valid_taxon])
    
    if (!is.null(taxa)) {
        taxa <- trimws(as.character(taxa))
        taxon2tips <- taxon2tips[names(taxon2tips) %in% taxa]
        
        if (length(taxon2tips) == 0L) {
            stop("None of the requested taxa were found.")
        }
    }
    
    if (length(taxon2tips) == 0L) {
        stop("No taxa remain to evaluate after filtering.")
    }
    
    # ---- 7) Internal helper: classify one taxon ------------------------------
    classify_one_taxon <- function(taxon_name, tip_idx) {
        n_tips <- length(tip_idx)
        
        if (n_tips == 1L) {
            return(data.frame(
                rank = rank_col,
                taxon = taxon_name,
                n_tips = n_tips,
                status = "singleton",
                largest_pure_cluster_n = 1L,
                largest_pure_cluster_frac = 1,
                second_pure_cluster_n = 0L,
                second_pure_cluster_frac = 0,
                n_pure_clusters = 1L,
                n_outside_largest_cluster = 0L,
                mrca_clade_size = 1L,
                mrca_contaminant_tips = 0L,
                mrca_purity = 1,
                interpretation = "Only one tip is present; monophyly is not informative.",
                stringsAsFactors = FALSE
            ))
        }
        
        if (n_tips < min_tips) {
            return(data.frame(
                rank = rank_col,
                taxon = taxon_name,
                n_tips = n_tips,
                status = "too_few_tips",
                largest_pure_cluster_n = NA_integer_,
                largest_pure_cluster_frac = NA_real_,
                second_pure_cluster_n = NA_integer_,
                second_pure_cluster_frac = NA_real_,
                n_pure_clusters = NA_integer_,
                n_outside_largest_cluster = NA_integer_,
                mrca_clade_size = NA_integer_,
                mrca_contaminant_tips = NA_integer_,
                mrca_purity = NA_real_,
                interpretation = "Too few tips to assess monophyly robustly.",
                stringsAsFactors = FALSE
            ))
        }
        
        target_count <- integer(total_nodes)
        target_count[tip_idx] <- 1L
        
        for (i in seq_len(nrow(edge))) {
            p <- edge[i, 1]
            ch <- edge[i, 2]
            target_count[p] <- target_count[p] + target_count[ch]
        }
        
        pure <- target_count > 0L & target_count == desc_total
        
        parent_is_pure <- rep(FALSE, total_nodes)
        has_parent <- !is.na(parent)
        parent_is_pure[has_parent] <- pure[parent[has_parent]]
        maximal_pure_nodes <- which(pure & !parent_is_pure)
        
        cluster_sizes <- sort(target_count[maximal_pure_nodes], decreasing = TRUE)
        n_pure_clusters <- length(cluster_sizes)
        
        largest_n <- if (n_pure_clusters >= 1L) cluster_sizes[1] else 0L
        second_n  <- if (n_pure_clusters >= 2L) cluster_sizes[2] else 0L
        
        largest_frac <- largest_n / n_tips
        second_frac <- second_n / n_tips
        outside_largest <- n_tips - largest_n
        
        mrca_node <- ape::getMRCA(tree2, tip_idx)
        if (is.null(mrca_node) || length(mrca_node) == 0L) {
            mrca_size <- n_tips
        } else {
            mrca_size <- desc_total[mrca_node]
        }
        
        mrca_contam <- mrca_size - n_tips
        mrca_purity <- n_tips / mrca_size
        
        monophyletic <- (n_pure_clusters == 1L && largest_n == n_tips)
        
        if (monophyletic) {
            status <- "monophyletic"
            interpretation <- "All members form one pure clade."
        } else if (largest_frac >= near_mono_frac &&
                   second_frac < polyphyly_second_frac) {
            status <- "near_monophyletic_minor_outliers"
            interpretation <- paste0(
                "Most members cluster together (largest pure cluster = ",
                round(100 * largest_frac, 1),
                "%), with only a few members outside; these may be mislabels, contaminants, ",
                "taxonomic outliers, or problematic placements."
            )
        } else if (largest_frac <= polyphyly_largest_max ||
                   second_frac >= polyphyly_second_frac) {
            status <- "polyphyletic"
            interpretation <- paste0(
                "Members are split into multiple substantial clusters (largest = ",
                round(100 * largest_frac, 1),
                "%, second = ",
                round(100 * second_frac, 1),
                "%), suggesting polyphyly at this rank."
            )
        } else {
            status <- "non_monophyletic_split"
            interpretation <- paste0(
                "The taxon is not monophyletic, but the split is intermediate: ",
                "not a single clean clade, and not strongly polyphyletic either."
            )
        }
        
        data.frame(
            rank = rank_col,
            taxon = taxon_name,
            n_tips = n_tips,
            status = status,
            largest_pure_cluster_n = largest_n,
            largest_pure_cluster_frac = largest_frac,
            second_pure_cluster_n = second_n,
            second_pure_cluster_frac = second_frac,
            n_pure_clusters = n_pure_clusters,
            n_outside_largest_cluster = outside_largest,
            mrca_clade_size = mrca_size,
            mrca_contaminant_tips = mrca_contam,
            mrca_purity = mrca_purity,
            interpretation = interpretation,
            stringsAsFactors = FALSE
        )
    }
    
    # ---- 8) Evaluate all taxa -----------------------------------------------
    taxon_names <- names(taxon2tips)
    
    if (verbose) {
        message("Evaluating ", length(taxon_names), " taxa at rank `", rank_col, "` ...")
    }
    
    out_list <- lapply(seq_along(taxon_names), function(i) {
        if (verbose && (i %% 50L == 0L || i == length(taxon_names))) {
            message("  processed ", i, " / ", length(taxon_names), " taxa")
        }
        classify_one_taxon(taxon_names[i], taxon2tips[[i]])
    })
    
    out <- do.call(rbind, out_list)
    rownames(out) <- NULL
    
    # ---- 9) Put the most problematic groups first ---------------------------
    status_levels <- c(
        "polyphyletic",
        "non_monophyletic_split",
        "near_monophyletic_minor_outliers",
        "monophyletic",
        "too_few_tips",
        "singleton"
    )
    
    out$status <- factor(out$status, levels = status_levels)
    out <- out[order(out$status, -out$n_tips, out$taxon), , drop = FALSE]
    out$status <- as.character(out$status)
    
    # ---- 10) Return tibble for tidyverse-friendly downstream use ------------
    tibble::as_tibble(out)
}
