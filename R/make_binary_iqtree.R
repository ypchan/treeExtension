#' Rebuild a rooted binary IQ-TREE by splitting the outgroup branch equally
#'
#' @description
#' Read an IQ-TREE Newick tree as a character string, identify the specified
#' outgroup tip and its branch length, remove the original outgroup occurrence,
#' and rebuild the root as:
#'
#' `(<ingroup>:L/2,<outgroup>:L/2);`
#'
#' The function performs direct Newick string manipulation and does not parse
#' or rewrite the tree using `ape`. Existing internal-node support labels are
#' therefore retained unchanged.
#'
#' @param tre_path Character scalar. Path to the input IQ-TREE Newick file.
#' @param out_lbl Character scalar. Exact tip label of the outgroup.
#'
#' @return
#' Invisibly returns the path to the generated `.bi.treefile`.
#'
#' @details
#' The outgroup must be attached directly to the current root and must contain
#' an explicit branch length.
#'
#' Supported branch-length formats include:
#'
#' - `0.125`
#' - `1.1111501163`
#' - `2`
#' - `.25`
#' - `1e-05`
#' - `1.25E+03`
#'
#' @examples
#' \dontrun{
#' make_binary_iqtree(
#'     tre_path = "example.treefile",
#'     out_lbl = "Outgroup_species"
#' )
#'
#' make_binary_iqtree(
#'     tre_path = "Dothidotthia.treefile",
#'     out_lbl = "Dothidotthia_robiniae_MFLUCC_16__1175"
#' )
#' }
#'
#' @export
make_binary_iqtree <- function(tre_path, out_lbl) {
    
    tre <- paste(
        readLines(tre_path, warn = FALSE),
        collapse = ""
    )
    
    if (!stringr::str_detect(
        tre,
        stringr::fixed(out_lbl)
    )) {
        stop(out_lbl, " not in tre_file")
    }
    
    # Support:
    #   0.125
    #   1.1111501163
    #   2
    #   .25
    #   1e-05
    #   1.25E+03
    branch_pattern <- paste0(
        "((?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)",
        "(?:[eE][+-]?[0-9]+)?)"
    )
    
    m <- stringr::str_match(
        tre,
        paste0(
            "\\Q",
            out_lbl,
            "\\E:",
            branch_pattern
        )
    )
    
    bl_str <- m[, 2]
    
    if (is.na(bl_str)) {
        stop(
            "No matching branch length found: ",
            out_lbl,
            ":[numeric branch length]"
        )
    }
    
    bl <- as.numeric(bl_str)
    
    if (is.na(bl)) {
        stop(
            "Extracted branch length could not be converted to numeric"
        )
    }
    
    tre2 <- stringr::str_replace(
        tre,
        stringr::fixed(
            paste0(
                ",",
                out_lbl,
                ":",
                bl_str
            )
        ),
        ""
    )
    
    if (identical(tre2, tre)) {
        tre2 <- stringr::str_replace(
            tre2,
            stringr::fixed(
                paste0(
                    out_lbl,
                    ":",
                    bl_str,
                    ","
                )
            ),
            ""
        )
    }
    
    if (identical(tre2, tre)) {
        stop(
            "No removable target strings found: ",
            paste0(",", out_lbl, ":", bl_str),
            " or ",
            paste0(out_lbl, ":", bl_str, ",")
        )
    }
    
    ingroup <- stringr::str_remove(
        tre2,
        stringr::fixed(";")
    )
    
    bl_half <- bl / 2
    
    bin_newick <- sprintf(
        "(%s:%s,%s:%s);",
        ingroup,
        bl_half,
        out_lbl,
        bl_half
    )
    
    out_file <- file.path(
        dirname(tre_path),
        paste0(
            tools::file_path_sans_ext(
                basename(tre_path)
            ),
            ".bi.treefile"
        )
    )
    
    writeLines(
        bin_newick,
        out_file
    )
    
    message(
        "Outgroup branch length: ",
        bl_str,
        "\nHalf branch length: ",
        bl_half,
        "\nOutput written to: ",
        out_file
    )
    
    invisible(out_file)
}