test_that("make_binary_iqtree returns and writes a rebuilt Newick tree", {
    tf <- tempfile(fileext = ".treefile")
    writeLines(
        "((A:0.1,B:0.2)0.9:0.3,C:0.4,Melanconis_stilbostoma_CFCC_50475:0.1023136003);",
        tf
    )
    
    res <- make_binary_iqtree(
        tf,
        "Melanconis_stilbostoma_CFCC_50475",
        digits = 8,
        return = "newick",
        write = TRUE
    )
    
    out_file <- file.path(
        dirname(tf),
        paste0(tools::file_path_sans_ext(basename(tf)), ".bi.treefile")
    )
    
    expect_true(file.exists(out_file))
    expect_equal(readLines(out_file), res)
    
    # Half of 0.1023136003 is 0.05115680015, formatted to 8 decimal places.
    expect_match(res, "Melanconis_stilbostoma_CFCC_50475:0\\.0511568")
    expect_match(
        res,
        "^\\(.*:0\\.0511568,Melanconis_stilbostoma_CFCC_50475:0\\.0511568\\);\\s*$"
    )
    
    expect_false(grepl("Melanconis_stilbostoma_CFCC_50475:0\\.1023136003", res))
    expect_false(grepl(",\\)", res, perl = TRUE))
    expect_false(grepl(",,", res, fixed = TRUE))
})

test_that("make_binary_iqtree works when outgroup is the first root child", {
    tf <- tempfile(fileext = ".treefile")
    writeLines(
        "(Melanconis_stilbostoma_CFCC_50475:0.8,D:0.2,(E:0.1,F:0.1):0.1);",
        tf
    )
    
    res <- make_binary_iqtree(
        tf,
        "Melanconis_stilbostoma_CFCC_50475",
        digits = 6,
        return = "newick"
    )
    
    expect_match(res, "Melanconis_stilbostoma_CFCC_50475:0\\.4")
    expect_match(
        res,
        "^\\(.*:0\\.4,Melanconis_stilbostoma_CFCC_50475:0\\.4\\);\\s*$"
    )
    expect_false(grepl(",,", res, fixed = TRUE))
    expect_false(grepl(",\\)", res, perl = TRUE))
})

test_that("make_binary_iqtree can return a phylo object", {
    phy <- ape::read.tree(text = "(A:0.2,B:0.4,Outgroup:0.6);")
    
    res <- make_binary_iqtree(
        phy,
        "Outgroup",
        return = "phylo"
    )
    
    expect_s3_class(res, "phylo")
    
    outgroup_tip <- match("Outgroup", res$tip.label)
    outgroup_edge <- which(res$edge[, 2] == outgroup_tip)
    expect_equal(res$edge.length[outgroup_edge], 0.3)
})

test_that("make_binary_iqtree rejects nested outgroups", {
    tf <- tempfile(fileext = ".treefile")
    writeLines("((Outgroup:0.5,A:0.1):0.2,B:0.3);", tf)
    
    expect_error(
        make_binary_iqtree(tf, "Outgroup", return = "newick"),
        "attached directly to the current root"
    )
})
