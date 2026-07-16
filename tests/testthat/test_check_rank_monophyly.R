test_that("check_rank_monophyly classifies simple monophyletic groups", {
    tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
    metadata <- data.frame(
        genome_label = c("A", "B", "C", "D"),
        family = c("Fam1", "Fam1", "Fam2", "Fam2")
    )
    
    res <- check_rank_monophyly(
        tree = tree,
        data = metadata,
        rank_col = "family",
        tip_col = "genome_label",
        check_rooted = FALSE
    )
    
    expect_s3_class(res, "tbl_df")
    status_by_taxon <- setNames(res$status, res$taxon)
    expect_equal(status_by_taxon[["Fam1"]], "monophyletic")
    expect_equal(status_by_taxon[["Fam2"]], "monophyletic")
})

test_that("check_rank_monophyly reports missing metadata clearly", {
    tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
    metadata <- data.frame(
        genome_label = c("A", "B"),
        family = c("Fam1", "Fam1")
    )
    
    expect_error(
        check_rank_monophyly(
            tree = tree,
            data = metadata,
            rank_col = "family",
            tip_col = "genome_label",
            check_rooted = FALSE
        ),
        "missing in `data`"
    )
})
