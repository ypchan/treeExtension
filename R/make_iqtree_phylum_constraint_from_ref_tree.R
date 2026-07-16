# ============================================================
# Build IQ-TREE phylum-level constraint tree from reference tree
#
# Purpose:
#   1. Read query sequence lineage table:
#        seq_id    lineage
#
#   2. Extract phylum for each query sequence.
#
#   3. Read a reference phylogenetic tree.
#      The reference tree can be:
#        A) species/genus/family/order-level tips + ref_lineage.tsv
#        B) phylum-level tips directly
#
#   4. Inherit phylum-level topology from reference tree.
#
#   5. Replace each phylum tip by a polytomy of query seq_id.
#
#   6. Output IQ-TREE constraint tree:
#        constraint.tree
#
# Usage in IQ-TREE:
#   iqtree2 -s aln.fa -m MFP -g constraint.tree -B 1000 -T AUTO --prefix constrained
#
# Important:
#   The final constraint tree tips MUST be seq_id in your alignment.
#   Taxon names / phylum names are only used to build the backbone.
# ============================================================

library(tidyverse)
library(ape)
library(this.path)
setwd('D:\\BaiduSyncdisk\\SZU\\postdoc\\Pacearchaeales\\part1_phylogenetics\\manuscript_hifi_16s\\supplementary_files\\gtdb232_tree_constraints')

# -----------------------------
# 1. User parameters
# -----------------------------

query_lineage_file <- "query_lineage.tsv"
ref_tree_file      <- "ref.tree"

# If reference tree tips are species/genus/family/order names, provide this file.
# It must contain columns:
#   taxon    lineage
#
# If reference tree tips are already phylum names, set:
#   ref_lineage_file <- NA
ref_lineage_file <- "ref_lineage.tsv"
# ref_lineage_file <- NA

out_constraint_tree <- "constraint.tree"
out_phylum_map      <- "query_seqid_phylum.tsv"
out_ref_backbone    <- "reference_phylum_backbone.tree"

# If TRUE, stop when a query phylum is not found in the reference tree.
# If FALSE, missing phyla are added as unresolved clades at the root.
strict_ref_phylum <- FALSE

# -----------------------------
# 2. Helper functions
# -----------------------------

clean_rank_name <- function(x) {
    x <- as.character(x)
    x <- str_trim(x)
    x <- str_replace(x, "^[a-zA-Z]__", "")
    x <- str_replace_all(x, "\\s+", "_")
    x
}

extract_phylum <- function(lineage) {
    lineage <- as.character(lineage)
    
    out <- map_chr(lineage, function(x) {
        parts <- str_split(x, ";", simplify = FALSE)[[1]]
        parts <- str_trim(parts)
        
        # Priority 1: detect p__ prefix
        p_hit <- parts[str_detect(parts, "^p__")]
        if (length(p_hit) > 0) {
            return(clean_rank_name(p_hit[1]))
        }
        
        # Priority 2: use second field as phylum if no p__ prefix
        # Example: Fungi;Ascomycota;Sordariomycetes
        if (length(parts) >= 2) {
            return(clean_rank_name(parts[2]))
        }
        
        return(NA_character_)
    })
    
    out
}

quote_newick_label <- function(x) {
    x <- as.character(x)
    
    # Newick-safe labels can be unquoted.
    safe <- str_detect(x, "^[A-Za-z0-9_.|+-]+$")
    
    x2 <- x
    x2[!safe] <- paste0("'", str_replace_all(x2[!safe], "'", "''"), "'")
    x2
}

strip_tree_edge_length <- function(tree) {
    tree$edge.length <- NULL
    tree
}

# Convert ape tree to Newick body, replacing tip labels by custom strings.
# replacement is a named vector/list:
#   names = phylum labels in backbone tree
#   values = sequence clades, e.g. "(seq1,seq2,seq3)"
tree_to_newick_body <- function(tree, replacement = NULL) {
    n_tip <- length(tree$tip.label)
    root_node <- n_tip + 1
    
    children <- split(tree$edge[, 2], tree$edge[, 1])
    
    rec <- function(node) {
        if (node <= n_tip) {
            lab <- tree$tip.label[node]
            
            if (!is.null(replacement) && lab %in% names(replacement)) {
                return(replacement[[lab]])
            } else {
                return(quote_newick_label(lab))
            }
        }
        
        ch <- children[[as.character(node)]]
        paste0("(", paste(map_chr(ch, rec), collapse = ","), ")")
    }
    
    rec(root_node)
}

make_seq_clade <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & ids != ""]
    ids <- sort(ids)
    
    if (length(ids) == 0) {
        return(NA_character_)
    }
    
    if (length(ids) == 1) {
        return(quote_newick_label(ids))
    }
    
    paste0("(", paste(quote_newick_label(ids), collapse = ","), ")")
}

normalize_tree_tip_labels <- function(tree) {
    tree$tip.label <- clean_rank_name(tree$tip.label)
    tree
}

# -----------------------------
# 3. Read query sequence lineage
# -----------------------------

query_tbl <- readr::read_tsv(query_lineage_file, show_col_types = FALSE)

if (!all(c("seq_id", "lineage") %in% colnames(query_tbl))) {
    stop("query_lineage.tsv must contain columns: seq_id, lineage")
}

query_tbl <- query_tbl %>%
    transmute(
        seq_id = as.character(seq_id),
        lineage = as.character(lineage),
        phylum = extract_phylum(lineage)
    ) %>%
    filter(!is.na(seq_id), seq_id != "") %>%
    filter(!is.na(phylum), phylum != "", phylum != "NA")

if (nrow(query_tbl) == 0) {
    stop("No valid query sequences with phylum found.")
}

readr::write_tsv(query_tbl, out_phylum_map)

seq_by_phylum <- split(query_tbl$seq_id, query_tbl$phylum)

query_phyla <- sort(names(seq_by_phylum))

message("[INFO] Query sequences: ", nrow(query_tbl))
message("[INFO] Query phyla: ", paste(query_phyla, collapse = ", "))

# Sequence clade for each query phylum
phylum_to_seq_clade <- map_chr(seq_by_phylum, make_seq_clade)

# -----------------------------
# 4. Read reference tree
# -----------------------------

ref_tree <- read.tree(ref_tree_file)
ref_tree <- strip_tree_edge_length(ref_tree)

if (is.na(ref_lineage_file)) {
    # ----------------------------------------------------------
    # Case A:
    # Reference tree tips are already phylum names.
    # ----------------------------------------------------------
    
    message("[INFO] Reference tree is treated as phylum-level tree.")
    
    ref_tree <- normalize_tree_tip_labels(ref_tree)
    
    ref_phyla <- sort(ref_tree$tip.label)
    
    common_phyla <- intersect(query_phyla, ref_phyla)
    missing_phyla <- setdiff(query_phyla, ref_phyla)
    
    if (length(common_phyla) == 0) {
        stop("No shared phyla between query data and reference phylum tree.")
    }
    
    if (length(missing_phyla) > 0) {
        msg <- paste0(
            "These query phyla are missing from reference tree: ",
            paste(missing_phyla, collapse = ", ")
        )
        
        if (strict_ref_phylum) {
            stop(msg)
        } else {
            warning(msg, "\nThey will be added as unresolved clades at the root.")
        }
    }
    
    phylum_tree <- keep.tip(ref_tree, common_phyla)
    phylum_tree <- strip_tree_edge_length(phylum_tree)
    
} else {
    # ----------------------------------------------------------
    # Case B:
    # Reference tree tips are taxon names.
    # Use ref_lineage.tsv to map reference tips to phylum.
    # ----------------------------------------------------------
    
    message("[INFO] Reference tree is treated as taxon-level tree.")
    message("[INFO] Reference lineage file: ", ref_lineage_file)
    
    ref_tax <- readr::read_tsv(ref_lineage_file, show_col_types = FALSE)
    
    if (!all(c("taxon", "lineage") %in% colnames(ref_tax))) {
        stop("ref_lineage.tsv must contain columns: taxon, lineage")
    }
    
    ref_tax <- ref_tax %>%
        transmute(
            taxon_raw = as.character(taxon),
            taxon = clean_rank_name(taxon),
            lineage = as.character(lineage),
            phylum = extract_phylum(lineage)
        ) %>%
        filter(!is.na(taxon), taxon != "") %>%
        filter(!is.na(phylum), phylum != "", phylum != "NA")
    
    ref_tree <- normalize_tree_tip_labels(ref_tree)
    
    ref_tax <- ref_tax %>%
        filter(taxon %in% ref_tree$tip.label)
    
    if (nrow(ref_tax) == 0) {
        stop("No ref_lineage.tsv taxon names match ref.tree tip labels.")
    }
    
    common_phyla <- intersect(query_phyla, unique(ref_tax$phylum))
    missing_phyla <- setdiff(query_phyla, unique(ref_tax$phylum))
    
    if (length(common_phyla) == 0) {
        stop("No shared phyla between query data and reference taxonomy.")
    }
    
    if (length(missing_phyla) > 0) {
        msg <- paste0(
            "These query phyla are missing from reference taxonomy/tree: ",
            paste(missing_phyla, collapse = ", ")
        )
        
        if (strict_ref_phylum) {
            stop(msg)
        } else {
            warning(msg, "\nThey will be added as unresolved clades at the root.")
        }
    }
    
    # Check phylum monophyly in reference tree.
    mono_check <- map_dfr(common_phyla, function(p) {
        tips <- ref_tax %>%
            filter(phylum == p) %>%
            pull(taxon) %>%
            unique()
        
        tips <- intersect(tips, ref_tree$tip.label)
        
        tibble(
            phylum = p,
            n_ref_tips = length(tips),
            is_monophyletic = if_else(
                length(tips) >= 2,
                is.monophyletic(ref_tree, tips),
                NA
            )
        )
    })
    
    readr::write_tsv(mono_check, "reference_phylum_monophyly_check.tsv")
    
    bad_mono <- mono_check %>%
        filter(n_ref_tips >= 2, is_monophyletic == FALSE)
    
    if (nrow(bad_mono) > 0) {
        warning(
            "Some phyla are not monophyletic in the reference tree: ",
            paste(bad_mono$phylum, collapse = ", "),
            "\nThe script will still choose one representative tip per phylum to inherit backbone topology."
        )
    }
    
    # Use one representative tip per phylum.
    # If each phylum is monophyletic, any representative preserves phylum-level backbone.
    rep_tbl <- ref_tax %>%
        filter(phylum %in% common_phyla) %>%
        group_by(phylum) %>%
        summarise(
            representative_tip = sort(unique(taxon))[1],
            n_ref_tips = n_distinct(taxon),
            .groups = "drop"
        )
    
    rep_tips <- rep_tbl$representative_tip
    
    phylum_tree <- keep.tip(ref_tree, rep_tips)
    phylum_tree <- strip_tree_edge_length(phylum_tree)
    
    # Rename representative species/taxon tips to phylum names.
    tip_map <- setNames(rep_tbl$phylum, rep_tbl$representative_tip)
    phylum_tree$tip.label <- tip_map[phylum_tree$tip.label]
}

# -----------------------------
# 5. Save reference phylum backbone
# -----------------------------

write.tree(phylum_tree, file = out_ref_backbone)
message("[DONE] ", out_ref_backbone)

# -----------------------------
# 6. Replace phylum tips by query seq_id clades
# -----------------------------

backbone_phyla <- phylum_tree$tip.label

replacement <- phylum_to_seq_clade[backbone_phyla]
names(replacement) <- backbone_phyla

constraint_body <- tree_to_newick_body(phylum_tree, replacement = replacement)

# Add missing query phyla as unresolved clades at the root.
missing_in_backbone <- setdiff(query_phyla, backbone_phyla)

if (length(missing_in_backbone) > 0) {
    missing_clades <- phylum_to_seq_clade[missing_in_backbone]
    constraint_body <- paste0(
        "(",
        constraint_body,
        ",",
        paste(missing_clades, collapse = ","),
        ")"
    )
}

constraint_tree <- paste0(constraint_body, ";")

writeLines(constraint_tree, out_constraint_tree)

message("[DONE] ", out_constraint_tree)

# -----------------------------
# 7. Report
# -----------------------------

cat("\n========== Summary ==========\n")
cat("Query sequences: ", nrow(query_tbl), "\n", sep = "")
cat("Query phyla: ", length(query_phyla), "\n", sep = "")
cat("Backbone phyla: ", length(backbone_phyla), "\n", sep = "")
cat("Missing phyla added at root: ", length(missing_in_backbone), "\n", sep = "")

if (length(missing_in_backbone) > 0) {
    cat("Missing phyla: ", paste(missing_in_backbone, collapse = ", "), "\n", sep = "")
}

cat("\nUse in IQ-TREE:\n")
cat("iqtree2 -s aln.fa -m MFP -g ", out_constraint_tree, " -B 1000 -T AUTO --prefix constrained\n", sep = "")