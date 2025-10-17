test_that("make_binary_iqtree splits outgroup length and roots binary", {
  tf <- tempfile(fileext = ".treefile")
  writeLines(
    "((A:0.1,B:0.2)0.9:0.3,C:0.4,Melanconis_stilbostoma_CFCC_50475:0.1023136003);",
    tf
  )

  out <- make_binary_iqtree(tf, "Melanconis_stilbostoma_CFCC_50475", digits = 8)
  expect_true(file.exists(out))

  res <- readLines(out)
  # Half of 0.1023136003 is 0.05115680015 -> formatted to 8 sig figs: 0.0511568
  expect_match(res, "^\\(.*Melanconis_stilbostoma_CFCC_50475:0\\.0511568\\):0\\.0511568;\\s*$")

  # Original token should be gone
  expect_false(grepl("\\(Melanconis_stilbostoma_CFCC_50475\\):0\\.1023136003", res, perl = TRUE))

  # Should not contain stray comma artifacts
  expect_false(grepl(",\\)", res, fixed = TRUE))
})

test_that("removal also works if outgroup is not last sibling (comma on right)", {
  tf <- tempfile(fileext = ".treefile")
  writeLines(
    "(Melanconis_stilbostoma_CFCC_50475:0.8,D:0.2,(E:0.1,F:0.1):0.1);",
    tf
  )

  out <- make_binary_iqtree(tf, "Melanconis_stilbostoma_CFCC_50475", digits = 6)
  res <- readLines(out)

  # Half of 0.8 is 0.4
  expect_match(res, ":0\\.4\\):0\\.4;\\s*$")
  expect_false(grepl(",,", res, fixed = TRUE))
})
