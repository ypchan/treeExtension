#' Map BLASTN (m6) pairwise similarities onto internal nodes of an IQ-TREE tree
#'
#' @description
#' Given a treedata object from `treeio::read.iqtree()` and a BLASTN outfmt 6 (m6)
#' file, this function enumerates all internal nodes, collects their descendant taxa,
#' and summarizes the **maximum** and **minimum** pairwise identities (pident)
#' within each node based on the best HSP per *undirected* taxa pair
#' (priority: highest `pident` -> longest `length` -> highest `bitscore`).
#'
#' It also prints diagnostics comparing taxa present in the tree vs. in the BLAST file.
#'
#' @param tre A treedata object returned by `treeio::read.iqtree("treefile")`.
#' @param blast_m6_path Path to a BLASTN outfmt 6 file (12-column tab-delimited).
#' @param len_threshold Numeric; minimum alignment `length` to retain (default 400).
#'
#' @return A tibble with one row per internal node containing:
#' \itemize{
#'   \item `node`: internal node id
#'   \item `ntaxa`: number of descendant taxa (that appear in BLAST)
#'   \item `taxa`: descendant taxa labels (semicolon-separated)
#'   \item `max_pident`, `max_pair_q`, `max_pair_s`, `max_length`, `max_gapopen`, `max_bitscore`
#'   \item `min_pident`, `min_pair_q`, `min_pair_s`, `min_length`, `min_gapopen`, `min_bitscore`
#' }
#'
#' @details
#' Self-hits (`qseqid == sseqid`) are removed. For each undirected taxa pair
#' (`A||B` equals `B||A`) only the **best** HSP is kept according to
#' `pident` (desc), then `length` (desc), then `bitscore` (desc).
#' Tip labels in the tree must match `qseqid/sseqid` in the m6 file.
#'
#' @examples
#' \dontrun{
#' tre <- treeio::read.iqtree("example.treefile")
#' res <- map_blastn_2_iqtree(
#'   tre,
#'   blast_m6_path = "blast.outfmt6",
#'   len_threshold = 400
#' )
#' head(res)
#' }
#'
#' @importFrom ggtree fortify
#' @importFrom tidytree offspring
#' @importFrom readr read_tsv
#' @importFrom dplyr filter pull transmute mutate across all_of left_join bind_rows arrange
#' @importFrom tibble tibble
#' @export
map_blastn_2_iqtree <- function(tre, blast_m6_path, len_threshold = 400) {
  # ---- 1) Tree to tibble via fortify() (keep this as requested) ----
  # Produces a tbl_tree-like data.frame with columns: node, parent, isTip, label, etc.
  tree_tbl <- ggtree::fortify(tre)
  internal_nodes <- dplyr::filter(tree_tbl, !isTip) |> dplyr::pull(node)
  tre_tips      <- dplyr::filter(tree_tbl,  isTip) |> dplyr::pull(label)

  # ---- 2) Read BLAST m6 and coerce numeric columns ----
  num_cols <- c("pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore")

  m6 <- readr::read_tsv(blast_m6_path, col_names = FALSE, show_col_types = FALSE) |>
    dplyr::transmute(
      qseqid   = X1,  sseqid   = X2,  pident  = X3,  length = X4,
      mismatch = X5,  gapopen  = X6,  qstart  = X7,  qend   = X8,
      sstart   = X9,  send     = X10, evalue  = X11, bitscore = X12
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(num_cols),
        ~ suppressWarnings(as.numeric(trimws(as.character(.))))
      )
    ) |>
    dplyr::filter(length >= len_threshold, qseqid != sseqid)

  # ---- 2.1) Diagnostics: taxa coverage ----
  m6_tips <- unique(c(m6$qseqid, m6$sseqid))
  diff12 <- setdiff(tre_tips, m6_tips)
  cat("Number of elements in tree: ", length(tre_tips), "\n")
  cat("Number of elements in blastn: ", length(m6_tips), "\n\n")
  cat("Set difference: ", paste(diff12, collapse = ", "), "\n")
  cat("Number of elements in this difference: ", length(diff12), "\n\n")

  # ---- 3) Undirected pairs; keep best HSP per pair (base R) ----
  # Best = highest pident -> longest alignment -> highest bitscore
  t1 <- pmin(m6$qseqid, m6$sseqid)                # undirected IDs (sorted)
  t2 <- pmax(m6$qseqid, m6$sseqid)
  pair_key <- paste(t1, t2, sep = "||")

  ord  <- order(-m6$pident, -m6$length, -m6$bitscore, na.last = TRUE)
  keep <- !duplicated(pair_key[ord])

  m6_pairs <- cbind(m6, t1 = t1, t2 = t2, pair_key = pair_key)[ord, ][keep, ]

  # Compact lookup for joins
  pair_lookup <- m6_pairs[, c("pair_key","pident","length","mismatch","gapopen",
                              "qseqid","sseqid","bitscore","qstart","qend","sstart","send","evalue")]

  # ---- 4) Helper: get descendant tips for a node (using tidytree::offspring) ----
  get_descendant_tips <- function(n) {
    desc <- tidytree::offspring(tree_tbl, n)
    if (nrow(desc) == 0) return(character(0))
    dplyr::filter(desc, isTip) |> dplyr::pull(label)
  }

  taxa_in_m6 <- unique(c(m6_pairs$t1, m6_pairs$t2))

  # ---- 5) Per-node min/max pident and associated stats ----
  res_list <- lapply(internal_nodes, function(n) {
    taxa <- get_descendant_tips(n)
    taxa <- intersect(taxa, taxa_in_m6)

    if (length(taxa) < 2) {
      return(tibble::tibble(
        node = n, ntaxa = length(taxa), taxa = paste(taxa, collapse = ";"),
        max_pident = NA_real_, max_pair_q = NA_character_, max_pair_s = NA_character_,
        max_length = NA_real_,  max_gapopen = NA_real_,  max_bitscore = NA_real_,
        min_pident = NA_real_, min_pair_q = NA_character_, min_pair_s = NA_character_,
        min_length = NA_real_,  min_gapopen = NA_real_,  min_bitscore = NA_real_
      ))
    }

    # enumerate all undirected within-node pairs
    combs <- t(combn(taxa, 2))
    pair_keys <- apply(combs, 1, function(x) {
      if (x[1] <= x[2]) paste(x[1], x[2], sep = "||") else paste(x[2], x[1], sep = "||")
    })

    df_pairs <- tibble::tibble(pair_key = pair_keys) |>
      dplyr::left_join(pair_lookup, by = "pair_key") |>
      dplyr::filter(!is.na(pident))

    if (nrow(df_pairs) == 0) {
      return(tibble::tibble(
        node = n, ntaxa = length(taxa), taxa = paste(taxa, collapse = ";"),
        max_pident = NA_real_, max_pair_q = NA_character_, max_pair_s = NA_character_,
        max_length = NA_real_,  max_gapopen = NA_real_,  max_bitscore = NA_real_,
        min_pident = NA_real_, min_pair_q = NA_character_, min_pair_s = NA_character_,
        min_length = NA_real_,  min_gapopen = NA_real_,  min_bitscore = NA_real_
      ))
    }

    i_max <- which.max(df_pairs$pident)
    i_min <- which.min(df_pairs$pident)
    row_max <- df_pairs[i_max, ]
    row_min <- df_pairs[i_min, ]

    tibble::tibble(
      node = n,
      ntaxa = length(taxa),
      taxa = paste(taxa, collapse = ";"),
      max_pident   = row_max$pident,
      max_pair_q   = row_max$qseqid,
      max_pair_s   = row_max$sseqid,
      max_length   = row_max$length,
      max_gapopen  = row_max$gapopen,
      max_bitscore = row_max$bitscore,
      min_pident   = row_min$pident,
      min_pair_q   = row_min$qseqid,
      min_pair_s   = row_min$sseqid,
      min_length   = row_min$length,
      min_gapopen  = row_min$gapopen,
      min_bitscore = row_min$bitscore
    )
  })

  dplyr::bind_rows(res_list) |> dplyr::arrange(node)
}
