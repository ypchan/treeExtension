make_binary_iqtree <- function(tre_path, out_lbl) {
  # Requires: stringr
  # install.packages("stringr") if needed
  library(stringr)

  # -- 1) Read and normalize ---------------------------------------------------
  # Read the Newick tree file and collapse into a single-line string.
  tre <- readLines(tre_path, warn = FALSE) |> paste(collapse = "")

  # Ensure the outgroup label exists (literal match, no regex semantics).
  if (!str_detect(tre, fixed(out_lbl))) stop(out_lbl, " not in tre_file")

  # -- 2) Parse branch length for outgroup ------------------------------------
  # Capture the branch length L immediately following "label:" (expects 0.xxx).
  # \Q...\E safely escapes any special characters in the label.
  m <- str_match(tre, paste0("\\Q", out_lbl, "\\E:(0\\.[0-9]+)"))
  bl_str <- m[, 2]                          # captured numeric string (e.g., "0.123")
  if (is.na(bl_str)) {
    stop("No matching pattern found: ", out_lbl, ":(0.[0-9]+)")
  }

  bl <- as.numeric(bl_str)                  # numeric branch length
  if (is.na(bl)) {
    stop("Extracted branch length could not be converted to numeric")
  }

  # -- 3) Remove the outgroup occurrence from the ingroup string --------------
  # Try removing ",label:L" first; if absent, try "label:L,".
  tre2 <- tre |> str_replace(fixed(paste0(",", out_lbl, ":", bl_str)), "")
  if (identical(tre2, tre)) {
    tre2 <- tre2 |> str_replace(fixed(paste0(out_lbl, ":", bl_str, ",")), "")
  }
  if (identical(tre2, tre)) {
    stop(
      "No removable target strings found: ",
      paste0(",", out_lbl, ":", bl_str), " or ",
      paste0(out_lbl, ":", bl_str, ",")
    )
  }

  # -- 4) Build the new binary tree -------------------------------------------
  # Remove trailing semicolon from ingroup and split the outgroup length equally.
  ingroup <- str_remove(tre2, fixed(";"))
  bl_half <- bl / 2
  bin_newick <- sprintf("(%s:%s,%s:%s);", ingroup, bl_half, out_lbl, bl_half)

  # -- 5) Emit output ----------------------------------------------------------
  # Write "<input>.bi.treefile" alongside the input file.
  out_file <- file.path(
    dirname(tre_path),
    paste0(tools::file_path_sans_ext(basename(tre_path)), ".bi.treefile")
  )
  writeLines(bin_newick, out_file)
  message("Output written to: ", out_file)
}