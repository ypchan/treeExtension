#' Build an IQ-TREE phylum-level constraint from a reference tree
#'
#' @description
#' Reads query sequence lineages, extracts their phyla, and replaces the
#' corresponding phylum tips in a reference-tree backbone with query-sequence
#' clades. Reference-tree tips may either be phylum names or lower-rank taxon
#' names accompanied by a lineage table.
#'
#' All input files are read only when this function is called. Loading or
#' installing `treeExtension` never changes the working directory or reads
#' analysis files.
#'
#' @param query_lineage_file Path to a tab-separated file containing columns
#'   `seq_id` and `lineage`.
#' @param ref_tree_file Path to the reference tree in Newick format.
#' @param ref_lineage_file Either `NA` when reference-tree tips are phylum
#'   names, or a path to a tab-separated file containing columns `taxon` and
#'   `lineage`.
#' @param out_constraint_tree Output path for the IQ-TREE constraint tree.
#' @param out_phylum_map Output path for the query sequence-to-phylum table.
#' @param out_ref_backbone Output path for the reference phylum backbone.
#' @param out_monophyly_check Output path for the reference phylum monophyly
#'   report. This file is written only when `ref_lineage_file` is supplied.
#' @param strict_ref_phylum If `TRUE`, stop when a query phylum is absent from
#'   the reference data. If `FALSE`, add absent phyla as unresolved root
#'   clades.
#' @param verbose If `TRUE`, print progress and summary messages.
#'
#' @return Invisibly returns a list containing output paths, a one-row summary,
#'   the query phylum map, and (when applicable) the monophyly report.
#'
#' @examples
#' \dontrun{
#' make_iqtree_phylum_constraint_from_ref_tree(
#'     query_lineage_file = "query_lineage.tsv",
#'     ref_tree_file = "ref.tree",
#'     ref_lineage_file = "ref_lineage.tsv"
#' )
#'
#' # Use NA when the reference tree already has phylum-level tips.
#' make_iqtree_phylum_constraint_from_ref_tree(
#'     query_lineage_file = "query_lineage.tsv",
#'     ref_tree_file = "phylum.tree",
#'     ref_lineage_file = NA
#' )
#' }
#'
#' @export
make_iqtree_phylum_constraint_from_ref_tree <- function(
        query_lineage_file,
        ref_tree_file,
        ref_lineage_file = NA_character_,
        out_constraint_tree = "constraint.tree",
        out_phylum_map = "query_seqid_phylum.tsv",
        out_ref_backbone = "reference_phylum_backbone.tree",
        out_monophyly_check = "reference_phylum_monophyly_check.tsv",
        strict_ref_phylum = FALSE,
        verbose = TRUE
) {
    .tree_extension_assert_input_file(query_lineage_file, "query_lineage_file")
    .tree_extension_assert_input_file(ref_tree_file, "ref_tree_file")

    if (length(ref_lineage_file) != 1L ||
        (!is.na(ref_lineage_file) && !is.character(ref_lineage_file))) {
        stop("`ref_lineage_file` must be one file path or NA.", call. = FALSE)
    }
    if (!is.na(ref_lineage_file)) {
        .tree_extension_assert_input_file(ref_lineage_file, "ref_lineage_file")
    } else {
        ref_lineage_file <- NA_character_
    }

    output_paths <- list(
        out_constraint_tree = out_constraint_tree,
        out_phylum_map = out_phylum_map,
        out_ref_backbone = out_ref_backbone,
        out_monophyly_check = out_monophyly_check
    )
    invisible(lapply(names(output_paths), function(name) {
        .tree_extension_assert_output_path(output_paths[[name]], name)
    }))

    if (!is.logical(strict_ref_phylum) ||
        length(strict_ref_phylum) != 1L || is.na(strict_ref_phylum)) {
        stop("`strict_ref_phylum` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    inform <- function(...) {
        if (verbose) {
            message(...)
        }
    }

    query_raw <- readr::read_tsv(query_lineage_file, show_col_types = FALSE)
    required_query_columns <- c("seq_id", "lineage")
    missing_query_columns <- setdiff(required_query_columns, names(query_raw))
    if (length(missing_query_columns) > 0L) {
        stop(
            "`query_lineage_file` is missing required columns: ",
            paste(missing_query_columns, collapse = ", "),
            call. = FALSE
        )
    }

    query_tbl <- data.frame(
        seq_id = as.character(query_raw$seq_id),
        lineage = as.character(query_raw$lineage),
        stringsAsFactors = FALSE
    )
    query_tbl$phylum <- .tree_extension_extract_phylum(query_tbl$lineage)
    keep_query <- !is.na(query_tbl$seq_id) & query_tbl$seq_id != "" &
        !is.na(query_tbl$phylum) & query_tbl$phylum != "" &
        query_tbl$phylum != "NA"
    query_tbl <- query_tbl[keep_query, , drop = FALSE]

    if (nrow(query_tbl) == 0L) {
        stop("No valid query sequences with a phylum were found.", call. = FALSE)
    }

    duplicated_seq_ids <- unique(query_tbl$seq_id[duplicated(query_tbl$seq_id)])
    if (length(duplicated_seq_ids) > 0L) {
        conflicting <- vapply(
            duplicated_seq_ids,
            function(id) length(unique(query_tbl$phylum[query_tbl$seq_id == id])) > 1L,
            logical(1)
        )
        if (any(conflicting)) {
            stop(
                "Some `seq_id` values map to more than one phylum: ",
                paste(duplicated_seq_ids[conflicting], collapse = ", "),
                call. = FALSE
            )
        }
        query_tbl <- query_tbl[!duplicated(query_tbl$seq_id), , drop = FALSE]
    }

    readr::write_tsv(query_tbl, out_phylum_map)

    seq_by_phylum <- split(query_tbl$seq_id, query_tbl$phylum)
    query_phyla <- sort(names(seq_by_phylum))
    phylum_to_seq_clade <- vapply(
        seq_by_phylum,
        .tree_extension_make_seq_clade,
        character(1)
    )

    inform("[INFO] Query sequences: ", nrow(query_tbl))
    inform("[INFO] Query phyla: ", paste(query_phyla, collapse = ", "))

    ref_tree <- ape::read.tree(ref_tree_file)
    if (!inherits(ref_tree, "phylo")) {
        stop("`ref_tree_file` did not contain one valid phylogenetic tree.", call. = FALSE)
    }
    ref_tree <- .tree_extension_strip_edge_length(ref_tree)
    mono_check <- NULL

    if (is.na(ref_lineage_file)) {
        inform("[INFO] Reference tree is treated as a phylum-level tree.")
        ref_tree <- .tree_extension_normalize_tree_tips(ref_tree)
        ref_phyla <- sort(ref_tree$tip.label)
        common_phyla <- intersect(query_phyla, ref_phyla)
        missing_phyla <- setdiff(query_phyla, ref_phyla)

        if (length(common_phyla) == 0L) {
            stop(
                "No shared phyla were found between the query data and reference tree.",
                call. = FALSE
            )
        }
        .tree_extension_handle_missing_phyla(
            missing_phyla,
            strict_ref_phylum,
            "reference tree"
        )

        phylum_tree <- ape::keep.tip(ref_tree, common_phyla)
        phylum_tree <- .tree_extension_strip_edge_length(phylum_tree)
    } else {
        inform("[INFO] Reference tree is treated as a taxon-level tree.")
        inform("[INFO] Reference lineage file: ", ref_lineage_file)

        ref_raw <- readr::read_tsv(ref_lineage_file, show_col_types = FALSE)
        required_ref_columns <- c("taxon", "lineage")
        missing_ref_columns <- setdiff(required_ref_columns, names(ref_raw))
        if (length(missing_ref_columns) > 0L) {
            stop(
                "`ref_lineage_file` is missing required columns: ",
                paste(missing_ref_columns, collapse = ", "),
                call. = FALSE
            )
        }

        ref_tax <- data.frame(
            taxon_raw = as.character(ref_raw$taxon),
            taxon = .tree_extension_clean_rank_name(ref_raw$taxon),
            lineage = as.character(ref_raw$lineage),
            stringsAsFactors = FALSE
        )
        ref_tax$phylum <- .tree_extension_extract_phylum(ref_tax$lineage)
        keep_ref <- !is.na(ref_tax$taxon) & ref_tax$taxon != "" &
            !is.na(ref_tax$phylum) & ref_tax$phylum != "" &
            ref_tax$phylum != "NA"
        ref_tax <- ref_tax[keep_ref, , drop = FALSE]

        ref_tree <- .tree_extension_normalize_tree_tips(ref_tree)
        ref_tax <- ref_tax[ref_tax$taxon %in% ref_tree$tip.label, , drop = FALSE]
        if (nrow(ref_tax) == 0L) {
            stop(
                "No taxon names in `ref_lineage_file` match the reference-tree tips.",
                call. = FALSE
            )
        }

        common_phyla <- intersect(query_phyla, unique(ref_tax$phylum))
        missing_phyla <- setdiff(query_phyla, unique(ref_tax$phylum))
        if (length(common_phyla) == 0L) {
            stop(
                "No shared phyla were found between the query data and reference taxonomy.",
                call. = FALSE
            )
        }
        .tree_extension_handle_missing_phyla(
            missing_phyla,
            strict_ref_phylum,
            "reference taxonomy/tree"
        )

        mono_rows <- lapply(common_phyla, function(phylum) {
            tips <- unique(ref_tax$taxon[ref_tax$phylum == phylum])
            tips <- intersect(tips, ref_tree$tip.label)
            data.frame(
                phylum = phylum,
                n_ref_tips = length(tips),
                is_monophyletic = if (length(tips) >= 2L) {
                    ape::is.monophyletic(ref_tree, tips)
                } else {
                    NA
                },
                stringsAsFactors = FALSE
            )
        })
        mono_check <- do.call(rbind, mono_rows)
        readr::write_tsv(mono_check, out_monophyly_check)

        bad_mono <- !is.na(mono_check$is_monophyletic) &
            !mono_check$is_monophyletic
        if (any(bad_mono)) {
            warning(
                "Some phyla are not monophyletic in the reference tree: ",
                paste(mono_check$phylum[bad_mono], collapse = ", "),
                ". One representative tip per phylum will be used for the backbone.",
                call. = FALSE
            )
        }

        ref_by_phylum <- split(ref_tax$taxon, ref_tax$phylum)
        ref_by_phylum <- ref_by_phylum[common_phyla]
        representative_tip <- vapply(
            ref_by_phylum,
            function(tips) sort(unique(tips))[1L],
            character(1)
        )
        rep_tbl <- data.frame(
            phylum = names(representative_tip),
            representative_tip = unname(representative_tip),
            n_ref_tips = vapply(ref_by_phylum, function(x) length(unique(x)), integer(1)),
            stringsAsFactors = FALSE
        )

        phylum_tree <- ape::keep.tip(ref_tree, rep_tbl$representative_tip)
        phylum_tree <- .tree_extension_strip_edge_length(phylum_tree)
        tip_map <- stats::setNames(rep_tbl$phylum, rep_tbl$representative_tip)
        phylum_tree$tip.label <- unname(tip_map[phylum_tree$tip.label])
    }

    ape::write.tree(phylum_tree, file = out_ref_backbone)
    inform("[DONE] ", out_ref_backbone)

    backbone_phyla <- phylum_tree$tip.label
    replacement <- phylum_to_seq_clade[backbone_phyla]
    names(replacement) <- backbone_phyla
    constraint_body <- .tree_extension_tree_to_newick_body(
        phylum_tree,
        replacement = replacement
    )

    missing_in_backbone <- setdiff(query_phyla, backbone_phyla)
    if (length(missing_in_backbone) > 0L) {
        missing_clades <- phylum_to_seq_clade[missing_in_backbone]
        constraint_body <- paste0(
            "(", constraint_body, ",",
            paste(missing_clades, collapse = ","), ")"
        )
    }

    constraint_tree <- paste0(constraint_body, ";")
    writeLines(constraint_tree, out_constraint_tree)
    inform("[DONE] ", out_constraint_tree)

    summary <- data.frame(
        query_sequences = nrow(query_tbl),
        query_phyla = length(query_phyla),
        backbone_phyla = length(backbone_phyla),
        missing_phyla_added_at_root = length(missing_in_backbone),
        stringsAsFactors = FALSE
    )

    if (verbose) {
        message(
            "[SUMMARY] Query sequences: ", summary$query_sequences,
            "; query phyla: ", summary$query_phyla,
            "; backbone phyla: ", summary$backbone_phyla,
            "; missing phyla added at root: ",
            summary$missing_phyla_added_at_root
        )
    }

    invisible(list(
        paths = c(
            constraint_tree = out_constraint_tree,
            phylum_map = out_phylum_map,
            reference_backbone = out_ref_backbone,
            monophyly_check = if (is.na(ref_lineage_file)) NA_character_ else out_monophyly_check
        ),
        summary = summary,
        query_phylum_map = tibble::as_tibble(query_tbl),
        monophyly_check = if (is.null(mono_check)) NULL else tibble::as_tibble(mono_check)
    ))
}

.tree_extension_assert_input_file <- function(path, argument) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
        stop("`", argument, "` must be one non-empty file path.", call. = FALSE)
    }
    if (!file.exists(path)) {
        stop("File supplied to `", argument, "` does not exist: ", path, call. = FALSE)
    }
}

.tree_extension_assert_output_path <- function(path, argument) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
        stop("`", argument, "` must be one non-empty file path.", call. = FALSE)
    }
}

.tree_extension_clean_rank_name <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^[A-Za-z]__", "", x)
    gsub("\\s+", "_", x)
}

.tree_extension_extract_phylum <- function(lineage) {
    vapply(as.character(lineage), function(value) {
        if (is.na(value)) {
            return(NA_character_)
        }
        parts <- trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
        phylum_part <- parts[grepl("^p__", parts)]
        if (length(phylum_part) > 0L) {
            return(.tree_extension_clean_rank_name(phylum_part[1L]))
        }
        if (length(parts) >= 2L) {
            return(.tree_extension_clean_rank_name(parts[2L]))
        }
        NA_character_
    }, character(1), USE.NAMES = FALSE)
}

.tree_extension_quote_newick_label <- function(x) {
    x <- as.character(x)
    safe <- grepl("^[A-Za-z0-9_.|+-]+$", x)
    x[!safe] <- paste0("'", gsub("'", "''", x[!safe], fixed = TRUE), "'")
    x
}

.tree_extension_strip_edge_length <- function(tree) {
    tree$edge.length <- NULL
    tree
}

.tree_extension_normalize_tree_tips <- function(tree) {
    tree$tip.label <- .tree_extension_clean_rank_name(tree$tip.label)
    tree
}

.tree_extension_make_seq_clade <- function(ids) {
    ids <- sort(unique(as.character(ids)))
    ids <- ids[!is.na(ids) & ids != ""]
    if (length(ids) == 0L) {
        return(NA_character_)
    }
    labels <- .tree_extension_quote_newick_label(ids)
    if (length(labels) == 1L) {
        return(labels)
    }
    paste0("(", paste(labels, collapse = ","), ")")
}

.tree_extension_tree_to_newick_body <- function(tree, replacement = NULL) {
    n_tip <- length(tree$tip.label)
    root_node <- n_tip + 1L
    children <- split(tree$edge[, 2L], tree$edge[, 1L])

    recurse <- function(node) {
        if (node <= n_tip) {
            label <- tree$tip.label[node]
            if (!is.null(replacement) && label %in% names(replacement)) {
                return(replacement[[label]])
            }
            return(.tree_extension_quote_newick_label(label))
        }
        child_nodes <- children[[as.character(node)]]
        paste0(
            "(",
            paste(vapply(child_nodes, recurse, character(1)), collapse = ","),
            ")"
        )
    }

    recurse(root_node)
}

.tree_extension_handle_missing_phyla <- function(
        missing_phyla,
        strict_ref_phylum,
        reference_description
) {
    if (length(missing_phyla) == 0L) {
        return(invisible(NULL))
    }
    message_text <- paste0(
        "These query phyla are missing from the ", reference_description, ": ",
        paste(missing_phyla, collapse = ", ")
    )
    if (strict_ref_phylum) {
        stop(message_text, call. = FALSE)
    }
    warning(
        message_text,
        ". They will be added as unresolved clades at the root.",
        call. = FALSE
    )
}
