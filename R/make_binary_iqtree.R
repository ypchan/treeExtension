#' Make a rooted binary IQ-TREE by splitting the outgroup branch equally
#'
#' @description
#' Read a Newick tree file produced by IQ-TREE, locate an outgroup tip by its
#' label, extract its pendant branch length `L`, remove the original occurrence,
#' and output a new binary tree `(<ingroup>:L/2,<outgroup>:L/2);`.
#'
#' @param tre_path Character, path to the input Newick tree file.
#' @param out_lbl  Character, the exact outgroup tip label to split on.
#'
#' @return Invisibly returns the output file path; writes a
#'   `*.bi.treefile` next to the input file.
#'
#' @details
#' This function expects the outgroup label to be present exactly once and its
#' pendant branch length to be of the form `:0.xxx`. The removal tries both
#' `,label:L` and `label:L,`. If neither matches, it errors with a diagnostic.
#'
#' @examples
#' \dontrun{
#' make_binary_iqtree("example.treefile", "Outgroup_Sp")
#' }
#'
#' @importFrom stringr str_detect str_match str_replace str_remove fixed
#' @importFrom tools file_path_sans_ext
#' @export
make_binary_iqtree <- function(tre_path, out_lbl) {
    # -- 1) Read and normalize ---------------------------------------------------
    tre <- readLines(tre_path, warn = FALSE) |> paste(collapse = "")
    
    # Ensure the outgroup label exists (literal match, no regex semantics).
    if (!stringr::str_detect(tre, stringr::fixed(out_lbl))) {
        stop(out_lbl, " not in tre_file")
    }
    
    # -- 2) Parse branch length for outgroup ------------------------------------
    # Capture the branch length L immediately following "label:" (expects 0.xxx).
    m <- stringr::str_match(tre, paste0("\\Q", out_lbl, "\\E:(0\\.[0-9]+)"))
    bl_str <- m[, 2]                          # e.g., "0.123"
    if (is.na(bl_str)) {
        stop("No matching pattern found: ", out_lbl, ":(0.[0-9]+)")
    }
    
    bl <- as.numeric(bl_str)
    if (is.na(bl)) stop("Extracted branch length could not be converted to numeric")
    
    # -- 3) Remove the outgroup occurrence from the ingroup string --------------
    tre2 <- tre |> stringr::str_replace(stringr::fixed(paste0(",", out_lbl, ":", bl_str)), "")
    if (identical(tre2, tre)) {
        tre2 <- tre2 |> stringr::str_replace(stringr::fixed(paste0(out_lbl, ":", bl_str, ",")), "")
    }
    if (identical(tre2, tre)) {
        stop("No removable target strings found: ",
             paste0(",", out_lbl, ":", bl_str), " or ",
             paste0(out_lbl, ":", bl_str, ","))
    }
    
    # -- 4) Build the new binary tree -------------------------------------------
    ingroup <- stringr::str_remove(tre2, stringr::fixed(";"))
    bl_half <- bl / 2
    bin_newick <- sprintf("(%s:%s,%s:%s);", ingroup, bl_half, out_lbl, bl_half)
    
    # -- 5) Emit output ----------------------------------------------------------
    out_file <- file.path(
        dirname(tre_path),
        paste0(tools::file_path_sans_ext(basename(tre_path)), ".bi.treefile")
    )
    writeLines(bin_newick, out_file)
    message("Output written to: ", out_file)
    invisible(out_file)
}
