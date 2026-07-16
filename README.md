# treeExtension

`treeExtension` is a small R package for practical phylogenetic-tree utilities:
rebuilding a rooted binary tree from an outgroup branch, checking monophyly
against metadata, and summarizing BLASTN pairwise similarity across internal
tree nodes.

The package is intentionally narrow. Each exported function does one job and
returns an R object that can be inspected, filtered, or written by the caller.

## Installation

```r
install.packages("remotes")
remotes::install_github("ypchan/treeExtension")
```

Some dependencies are Bioconductor packages. If installation reports missing
`ggtree`, `treeio`, or `tidytree`, install them first:

```r
install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio", "tidytree"))
```

## Functions

| Function | Purpose | Main output |
| --- | --- | --- |
| `make_binary_iqtree()` | Rebuild a rooted binary tree by splitting a root-attached outgroup branch into `L/2 + L/2` | `treedata`, `phylo`, or Newick string |
| `check_rank_monophyly()` | Check whether taxa at a metadata rank form clean clades on a rooted tree | tibble, one row per taxon |
| `summarize_blastn_by_node()` | Summarize observed BLASTN row-level identity statistics for each internal node | tibble, one row per internal node |

## Rebuild a Binary Outgroup Root

`make_binary_iqtree()` takes a tree file path, a `phylo` object, or a `treedata`
object. The outgroup tip must already be attached directly to the current root;
the function does not search for or reroot by a nested outgroup.

```r
library(treeExtension)

td <- make_binary_iqtree(
    tree = "example.treefile",
    outgroup_label = "Outgroup_Sp"
)
```

Return a Newick string instead:

```r
nwk <- make_binary_iqtree(
    tree = "example.treefile",
    outgroup_label = "Outgroup_Sp",
    return = "newick",
    digits = 8
)
```

Write the rebuilt tree to disk:

```r
make_binary_iqtree(
    tree = "example.treefile",
    outgroup_label = "Outgroup_Sp",
    write = TRUE,
    file = "example.bi.treefile",
    overwrite = TRUE
)
```

Important assumptions:

- The outgroup label must match exactly one tip.
- The outgroup branch must have a finite, non-negative branch length.
- The outgroup must be a direct child of the current root.
- For `treedata` inputs, topology and branch lengths are rebuilt; auxiliary
  node or tip annotations are not reattached after the topology change.

## Check Rank-Level Monophyly

`check_rank_monophyly()` aligns metadata to `tree$tip.label` using `tip_col`,
then evaluates one metadata rank at a time.

```r
res_family <- check_rank_monophyly(
    tree = tre,
    data = metadata,
    rank_col = "family",
    tip_col = "genome_label"
)

table(res_family$status)
subset(res_family, status != "monophyletic")
```

Returned columns include:

- `taxon`, `n_tips`, and `status`
- largest and second-largest pure cluster sizes
- MRCA clade size, contaminant count, and purity
- a short `interpretation`

Status values:

- `monophyletic`: all members form one pure clade
- `near_monophyletic_minor_outliers`: one dominant clean cluster plus minor
  outliers
- `polyphyletic`: multiple substantial clusters
- `non_monophyletic_split`: intermediate non-monophyly
- `singleton` or `too_few_tips`: not enough tips for robust interpretation

## Summarize BLASTN by Internal Node

`summarize_blastn_by_node()` accepts a `phylo` or `treedata` tree and BLASTN
outfmt 6 data. `blast` can be a file path or a data.frame/tibble containing the
standard first 12 outfmt 6 columns:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

```r
tre <- treeio::read.iqtree("example.treefile")

node_blast <- summarize_blastn_by_node(
    tree = tre,
    blast = "blast.outfmt6",
    min_alignment_length = 400,
    pair_policy = "all_hsps",
    diagnostics = TRUE
)

head(node_blast)
attr(node_blast, "diagnostics")
```

By default, every filtered BLASTN row/HSP is considered when computing the
within-node minimum and maximum identity. For a given internal node, the
function uses BLASTN rows whose query and subject are both descendant taxa of
that node. The output records the minimum and maximum `pident` values and the
query/subject taxa where those values occur.

Output columns include `node`, `ntaxa`, `taxa`, `n_blast_rows`,
`n_unique_pairs`, `min_pident`, `min_pair_q`, `min_pair_s`, `max_pident`,
`max_pair_q`, and `max_pair_s`.

Use `pair_policy = "best_hsp"` if multiple HSPs for the same undirected pair
should first be collapsed to the best HSP using this priority:

1. highest `pident`
2. longest `length`
3. highest `bitscore`

The function summarizes observed BLAST pairs only. It removes self-hits and
rows with missing query ID, subject ID, alignment length, or percent identity.
It does not infer missing pairwise comparisons.

## Development

Run these commands from the package root:

```r
library(devtools)

devtools::document()
devtools::load_all()
devtools::test()
devtools::check()
```

Before exporting a new function:

- keep the function focused on tree, metadata, IQ-TREE, or node-level summary
  work
- validate inputs early with specific error messages
- return a stable R object
- add roxygen documentation
- add tests under `tests/testthat/`
- update this README if the user-facing API changes

## License

MIT
