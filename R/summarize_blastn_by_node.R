#' Summarize BLASTN pairwise similarity statistics for internal tree nodes
#'
#' @description
#' Given a phylogenetic tree and BLASTN outfmt 6 results, summarize observed
#' within-node pairwise identity statistics for descendant taxa at each internal
#' node.
#'
#' By default, every filtered BLASTN row/HSP whose query and subject are both
#' descendant taxa of a node is considered in that node's min/max summary. If
#' `pair_policy = "best_hsp"`, each undirected taxa pair is first collapsed to
#' its best HSP using the ranking: highest `pident`, then longest alignment,
#' then highest `bitscore`.
#'
#' @param tree A `treedata` or `phylo` object.
#' @param blast Either a path to a BLASTN outfmt 6 file, or a data.frame/tibble
#'   containing the 12 standard BLAST outfmt 6 columns.
#' @param min_alignment_length Numeric scalar. Minimum alignment length retained.
#'   Default: `400`.
#' @param pair_policy Character scalar. One of `"all_hsps"` or `"best_hsp"`.
#'   Default: `"all_hsps"`.
#' @param diagnostics Logical. If `TRUE`, attach diagnostic information as an
#'   attribute named `"diagnostics"`. Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print informative messages.
#'   Default: `FALSE`.
#'
#' @return
#' A tibble with one row per internal node and summary statistics of observed
#' BLASTN pairwise similarity among descendant taxa.
#'
#' The returned tibble contains these columns:
#' - `node`
#' - `ntaxa`: number of descendant taxa in the node
#' - `taxa`: descendant taxa labels
#' - `n_blast_rows`: number of filtered BLASTN rows used for the node
#' - `n_unique_pairs`: number of unique undirected taxa pairs represented by
#'   those rows
#' - `max_pident`, `max_pair_q`, `max_pair_s`, `max_length`, `max_gapopen`,
#'   `max_bitscore`
#' - `min_pident`, `min_pair_q`, `min_pair_s`, `min_length`, `min_gapopen`,
#'   `min_bitscore`
#'
#' @details
#' Self-hits and rows with missing query IDs, subject IDs, alignment length, or
#' percent identity are removed. The function summarizes observed BLAST pairs;
#' it does not impute missing pairwise comparisons.
#'
#' @examples
#' \dontrun{
#' tre <- treeio::read.iqtree("example.treefile")
#' res <- summarize_blastn_by_node(
#'   tree = tre,
#'   blast = "blast.outfmt6"
#' )
#' }
#'
#' @export
summarize_blastn_by_node <- function(
        tree,
        blast,
        min_alignment_length = 400,
        pair_policy = c("all_hsps", "best_hsp"),
        diagnostics = FALSE,
        verbose = FALSE
) {
    pair_policy <- match.arg(pair_policy)
    
    # ---- 1) Validate input ---------------------------------------------------
    if (!inherits(tree, "treedata") && !inherits(tree, "phylo")) {
        stop("`tree` must be a `treedata` or `phylo` object.")
    }
    
    if (!is.numeric(min_alignment_length) ||
        length(min_alignment_length) != 1L ||
        is.na(min_alignment_length) ||
        min_alignment_length < 0) {
        stop("`min_alignment_length` must be a single non-negative numeric value.")
    }
    
    if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
        stop("`diagnostics` must be TRUE or FALSE.")
    }
    
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.")
    }
    
    # ---- 2) Fortify tree and collect nodes ----------------------------------
    tree_tbl <- as.data.frame(ggtree::fortify(tree))
    
    required_tree_cols <- c("node", "isTip", "label")
    missing_tree_cols <- setdiff(required_tree_cols, colnames(tree_tbl))
    if (length(missing_tree_cols) > 0L) {
        stop(
            "Fortified tree is missing required columns: ",
            paste(missing_tree_cols, collapse = ", ")
        )
    }
    
    internal_nodes <- tree_tbl$node[!tree_tbl$isTip]
    tree_tips <- tree_tbl$label[tree_tbl$isTip]
    
    # ---- 3) Standardize BLAST input -----------------------------------------
    standardize_blast_m6 <- function(x) {
        std_cols <- c(
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        )
        num_cols <- setdiff(std_cols, c("qseqid", "sseqid"))
        
        if (is.character(x) && length(x) == 1L) {
            if (!file.exists(x)) {
                stop("BLAST file does not exist: ", x)
            }
            
            raw_m6 <- readr::read_tsv(
                x,
                col_names = FALSE,
                show_col_types = FALSE
            )
            
            if (ncol(raw_m6) < 12L) {
                stop("BLAST outfmt 6 file must contain at least 12 columns.")
            }
            
            raw_m6 <- as.data.frame(raw_m6)
            m6 <- raw_m6[, seq_len(12L), drop = FALSE]
            colnames(m6) <- std_cols
        } else if (is.data.frame(x)) {
            m6 <- as.data.frame(x)
            
            if (all(std_cols %in% colnames(m6))) {
                m6 <- m6[, std_cols, drop = FALSE]
            } else if (ncol(m6) >= 12L) {
                m6 <- m6[, seq_len(12L), drop = FALSE]
                colnames(m6) <- std_cols
            } else {
                stop(
                    "`blast` data.frame must either contain the standard BLAST outfmt 6 ",
                    "column names, or have at least 12 columns."
                )
            }
        } else {
            stop("`blast` must be either a file path or a data.frame/tibble.")
        }
        
        m6$qseqid <- as.character(m6$qseqid)
        m6$sseqid <- as.character(m6$sseqid)
        
        for (col in num_cols) {
            m6[[col]] <- suppressWarnings(as.numeric(trimws(as.character(m6[[col]]))))
        }
        
        keep <- !is.na(m6$qseqid) &
            !is.na(m6$sseqid) &
            m6$qseqid != "" &
            m6$sseqid != "" &
            m6$qseqid != m6$sseqid &
            !is.na(m6$length) &
            m6$length >= min_alignment_length &
            !is.na(m6$pident)
        
        tibble::as_tibble(m6[keep, , drop = FALSE])
    }
    
    m6 <- standardize_blast_m6(blast)
    
    # ---- 4) Coverage diagnostics --------------------------------------------
    blast_tips <- unique(c(m6$qseqid, m6$sseqid))
    missing_in_blast <- setdiff(tree_tips, blast_tips)
    missing_in_tree <- setdiff(blast_tips, tree_tips)
    
    if (verbose) {
        message("Number of tips in tree: ", length(tree_tips))
        message("Number of unique taxa in BLAST: ", length(blast_tips))
        message("Tips present in tree but absent from BLAST: ", length(missing_in_blast))
        message("Tips present in BLAST but absent from tree: ", length(missing_in_tree))
    }
    
    # ---- 5) Apply pair policy -----------------------------------------------
    m6$t1 <- pmin(m6$qseqid, m6$sseqid)
    m6$t2 <- pmax(m6$qseqid, m6$sseqid)
    m6$pair_key <- paste(m6$t1, m6$t2, sep = "||")
    
    if (identical(pair_policy, "best_hsp")) {
        ord <- order(-m6$pident, -m6$length, -m6$bitscore, na.last = TRUE)
        keep_best <- ord[!duplicated(m6$pair_key[ord])]
        m6_pairs <- m6[keep_best, , drop = FALSE]
    } else {
        m6_pairs <- m6
    }
    
    # ---- 6) Helper: get descendant tips of a node ---------------------------
    get_descendant_tips <- function(node_id) {
        desc <- as.data.frame(tidytree::offspring(tree_tbl, node_id))
        if (nrow(desc) == 0L) {
            return(character(0))
        }
        
        desc$label[desc$isTip]
    }
    
    empty_node_summary <- function(node_id, taxa) {
        tibble::tibble(
            node = node_id,
            ntaxa = length(taxa),
            taxa = paste(taxa, collapse = ";"),
            n_blast_rows = 0L,
            n_unique_pairs = 0L,
            max_pident = NA_real_,
            max_pair_q = NA_character_,
            max_pair_s = NA_character_,
            max_length = NA_real_,
            max_gapopen = NA_real_,
            max_bitscore = NA_real_,
            min_pident = NA_real_,
            min_pair_q = NA_character_,
            min_pair_s = NA_character_,
            min_length = NA_real_,
            min_gapopen = NA_real_,
            min_bitscore = NA_real_
        )
    }

    select_extreme_identity_row <- function(pair_rows, maximum) {
        if (maximum) {
            ord <- order(
                -pair_rows$pident,
                -pair_rows$length,
                -pair_rows$bitscore,
                pair_rows$t1,
                pair_rows$t2,
                na.last = TRUE
            )
        } else {
            ord <- order(
                pair_rows$pident,
                -pair_rows$length,
                -pair_rows$bitscore,
                pair_rows$t1,
                pair_rows$t2,
                na.last = TRUE
            )
        }

        pair_rows[ord[1L], , drop = FALSE]
    }
    
    # ---- 7) Per-node summary -------------------------------------------------
    summarize_one_node <- function(node_id) {
        taxa <- get_descendant_tips(node_id)
        
        if (length(taxa) < 2L) {
            return(empty_node_summary(node_id, taxa))
        }
        
        pair_rows <- m6_pairs[
            m6_pairs$t1 %in% taxa & m6_pairs$t2 %in% taxa,
            ,
            drop = FALSE
        ]
        
        if (nrow(pair_rows) == 0L) {
            return(empty_node_summary(node_id, taxa))
        }
        
        row_max <- select_extreme_identity_row(pair_rows, maximum = TRUE)
        row_min <- select_extreme_identity_row(pair_rows, maximum = FALSE)
        
        tibble::tibble(
            node = node_id,
            ntaxa = length(taxa),
            taxa = paste(taxa, collapse = ";"),
            n_blast_rows = nrow(pair_rows),
            n_unique_pairs = length(unique(pair_rows$pair_key)),
            max_pident = row_max$pident[[1]],
            max_pair_q = row_max$qseqid[[1]],
            max_pair_s = row_max$sseqid[[1]],
            max_length = row_max$length[[1]],
            max_gapopen = row_max$gapopen[[1]],
            max_bitscore = row_max$bitscore[[1]],
            min_pident = row_min$pident[[1]],
            min_pair_q = row_min$qseqid[[1]],
            min_pair_s = row_min$sseqid[[1]],
            min_length = row_min$length[[1]],
            min_gapopen = row_min$gapopen[[1]],
            min_bitscore = row_min$bitscore[[1]]
        )
    }
    
    if (verbose) {
        message("Summarizing BLASTN statistics for ", length(internal_nodes), " internal nodes ...")
    }
    
    out_list <- lapply(seq_along(internal_nodes), function(i) {
        if (verbose && (i %% 100L == 0L || i == length(internal_nodes))) {
            message("  processed ", i, " / ", length(internal_nodes), " nodes")
        }
        summarize_one_node(internal_nodes[i])
    })
    
    out <- tibble::as_tibble(do.call(rbind, out_list))
    out <- out[order(out$node), , drop = FALSE]
    
    # ---- 8) Attach diagnostics optionally -----------------------------------
    if (diagnostics) {
        attr(out, "diagnostics") <- list(
            n_tree_tips = length(tree_tips),
            n_blast_tips = length(blast_tips),
            missing_in_blast = missing_in_blast,
            missing_in_tree = missing_in_tree,
            pair_policy = pair_policy,
            min_alignment_length = min_alignment_length
        )
    }
    
    out
}
