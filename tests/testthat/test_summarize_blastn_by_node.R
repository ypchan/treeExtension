test_that("summarize_blastn_by_node summarizes all observed internal-node HSP rows by default", {
    tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
    blast <- data.frame(
        qseqid = c("A", "A", "A", "B"),
        sseqid = c("B", "B", "C", "C"),
        pident = c(98, 99, 90, 91),
        length = c(600, 500, 700, 650),
        mismatch = c(2, 1, 10, 9),
        gapopen = c(0, 0, 1, 1),
        qstart = c(1, 1, 1, 1),
        qend = c(600, 500, 700, 650),
        sstart = c(1, 1, 1, 1),
        send = c(600, 500, 700, 650),
        evalue = c(0, 0, 1e-20, 1e-20),
        bitscore = c(800, 900, 700, 710)
    )
    
    res <- summarize_blastn_by_node(
        tree = tree,
        blast = blast,
        min_alignment_length = 400
    )
    
    expect_s3_class(res, "tbl_df")
    taxa_key <- vapply(
        strsplit(res$taxa, ";", fixed = TRUE),
        function(x) paste(sort(x), collapse = ";"),
        character(1)
    )
    expect_true(any(taxa_key == "A;B"))
    
    ab_row <- res[taxa_key == "A;B", , drop = FALSE]
    expect_equal(ab_row$ntaxa, 2)
    expect_equal(ab_row$n_blast_rows, 2)
    expect_equal(ab_row$n_unique_pairs, 1)
    expect_equal(ab_row$max_pident, 99)
    expect_equal(ab_row$min_pident, 98)
    expect_equal(ab_row$max_pair_q, "A")
    expect_equal(ab_row$max_pair_s, "B")
    expect_equal(ab_row$min_pair_q, "A")
    expect_equal(ab_row$min_pair_s, "B")
    
    root_row <- res[taxa_key == "A;B;C", , drop = FALSE]
    expect_equal(root_row$n_blast_rows, 4)
    expect_equal(root_row$n_unique_pairs, 3)
    expect_equal(root_row$max_pident, 99)
    expect_equal(root_row$min_pident, 90)
    expect_equal(root_row$min_pair_q, "A")
    expect_equal(root_row$min_pair_s, "C")
})

test_that("summarize_blastn_by_node can collapse each pair to its best HSP", {
    tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
    blast <- data.frame(
        qseqid = c("A", "A", "A", "B"),
        sseqid = c("B", "B", "C", "C"),
        pident = c(98, 99, 90, 91),
        length = c(600, 500, 700, 650),
        mismatch = c(2, 1, 10, 9),
        gapopen = c(0, 0, 1, 1),
        qstart = c(1, 1, 1, 1),
        qend = c(600, 500, 700, 650),
        sstart = c(1, 1, 1, 1),
        send = c(600, 500, 700, 650),
        evalue = c(0, 0, 1e-20, 1e-20),
        bitscore = c(800, 900, 700, 710)
    )
    
    res <- summarize_blastn_by_node(
        tree = tree,
        blast = blast,
        min_alignment_length = 400,
        pair_policy = "best_hsp"
    )
    
    taxa_key <- vapply(
        strsplit(res$taxa, ";", fixed = TRUE),
        function(x) paste(sort(x), collapse = ";"),
        character(1)
    )
    
    ab_row <- res[taxa_key == "A;B", , drop = FALSE]
    expect_equal(ab_row$n_blast_rows, 1)
    expect_equal(ab_row$n_unique_pairs, 1)
    expect_equal(ab_row$max_pident, 99)
    expect_equal(ab_row$min_pident, 99)
    
    root_row <- res[taxa_key == "A;B;C", , drop = FALSE]
    expect_equal(root_row$n_blast_rows, 3)
    expect_equal(root_row$n_unique_pairs, 3)
    expect_equal(root_row$max_pident, 99)
    expect_equal(root_row$min_pident, 90)
})

test_that("summarize_blastn_by_node can attach diagnostics", {
    tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
    blast <- data.frame(
        qseqid = "A",
        sseqid = "B",
        pident = 99,
        length = 500,
        mismatch = 1,
        gapopen = 0,
        qstart = 1,
        qend = 500,
        sstart = 1,
        send = 500,
        evalue = 0,
        bitscore = 900
    )
    
    res <- summarize_blastn_by_node(
        tree = tree,
        blast = blast,
        diagnostics = TRUE
    )
    
    diag <- attr(res, "diagnostics")
    expect_equal(diag$n_tree_tips, 3)
    expect_equal(diag$missing_in_blast, "C")
    
    taxa_key <- vapply(
        strsplit(res$taxa, ";", fixed = TRUE),
        function(x) paste(sort(x), collapse = ";"),
        character(1)
    )
    root_row <- res[taxa_key == "A;B;C", , drop = FALSE]
    expect_equal(root_row$ntaxa, 3)
    expect_equal(root_row$n_blast_rows, 1)
})
