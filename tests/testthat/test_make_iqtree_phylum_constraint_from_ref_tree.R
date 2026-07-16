test_that("constraint inputs are read only when the function is called", {
    expect_true(is.function(make_iqtree_phylum_constraint_from_ref_tree))
})

test_that("a phylum-level reference tree produces a query-tip constraint", {
    skip_if_not_installed("ape")
    skip_if_not_installed("readr")

    work_dir <- tempfile("treeExtension-constraint-")
    dir.create(work_dir)

    query_file <- file.path(work_dir, "query.tsv")
    ref_tree_file <- file.path(work_dir, "ref.tree")
    constraint_file <- file.path(work_dir, "constraint.tree")
    map_file <- file.path(work_dir, "map.tsv")
    backbone_file <- file.path(work_dir, "backbone.tree")

    writeLines(
        c(
            "seq_id\tlineage",
            "seq 1\td__Bacteria;p__Firmicutes;c__Bacilli",
            "seq2\td__Bacteria;p__Firmicutes;c__Bacilli",
            "seq3\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria"
        ),
        query_file
    )
    writeLines("(p__Firmicutes,p__Proteobacteria);", ref_tree_file)

    result <- make_iqtree_phylum_constraint_from_ref_tree(
        query_lineage_file = query_file,
        ref_tree_file = ref_tree_file,
        ref_lineage_file = NA,
        out_constraint_tree = constraint_file,
        out_phylum_map = map_file,
        out_ref_backbone = backbone_file,
        out_monophyly_check = file.path(work_dir, "unused.tsv"),
        verbose = FALSE
    )

    expect_true(all(file.exists(constraint_file, map_file, backbone_file)))
    constraint <- paste(readLines(constraint_file), collapse = "")
    expect_match(constraint, "'seq 1'", fixed = TRUE)
    expect_match(constraint, "seq2", fixed = TRUE)
    expect_match(constraint, "seq3", fixed = TRUE)
    expect_false(grepl("Firmicutes|Proteobacteria", constraint))
    expect_equal(result$summary$query_sequences, 3)
    expect_equal(result$summary$query_phyla, 2)
})

test_that("file errors identify the function argument", {
    expect_error(
        make_iqtree_phylum_constraint_from_ref_tree(
            query_lineage_file = tempfile("missing-query-"),
            ref_tree_file = tempfile("missing-tree-")
        ),
        "query_lineage_file"
    )
})
