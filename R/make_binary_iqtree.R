#' Convert IQ-TREE .treefile into a rooted binary tree by splitting the outgroup branch (L/2 + L/2)
#'
#' @description
#' Text-only (base R) transformation supporting exactly two token forms:
#'   1) "(<label>:<number>"   — immediately after an opening parenthesis
#'   2) ",<label>:<number>"   — immediately after a comma, with **no spaces** after the comma
#'
#' Steps:
#'   - Find the first occurrence of the outgroup token in either of the above forms.
#'   - Extract the pendant branch length L.
#'   - Remove the token from the original Newick string while preserving comma balance
#'     (consume exactly one neighboring comma if present to avoid creating ",)" or ",,").
#'   - Wrap the remaining ingroup and the outgroup into a new outer root, assigning L/2
#'     to the outgroup edge and L/2 to the root-to-ingroup edge:
#'         ( <ingroup> , <outgroup_label>:L/2 ):L/2;
#'   - Write the result to "<input>.bi.treefile".
#'
#' @param iqtree_treefile character(1). Path to input IQ-TREE `.treefile` (single-line Newick).
#' @param outgroup_label  character(1). Exact tip label to be used as outgroup (no quotes).
#' @param digits          NULL or integer(1). If not NULL, print L/2 with this many significant
#'                        digits via `formatC(format = "g")` (stable, concise output).
#'
#' @return Invisibly returns the output path `<input>.bi.treefile`.
#' @export
#'
#' @examples
#' \dontrun{
#' make_binary_iqtree("05_iqtree/iqtree_ml.treefile", "Diaporthella_corylina_CBS_121124")
#' }
make_binary_iqtree <- function(iqtree_treefile, outgroup_label, digits = NULL) {
  ## --------------------------- Local helpers (base R only) ---------------------------

  # Escape regex metacharacters so the literal label is matched exactly.
  re_escape <- function(x) gsub("([][{}()+*^$|?\\.^-])", "\\\\\\\\1", x, perl = TRUE)

  # Read the entire file as a single string (IQ-TREE usually writes one-line Newick).
  read_one_line <- function(path) paste(readLines(path, warn = FALSE), collapse = "")

  # Format L/2 either as raw numeric (default) or with controlled significant digits.
  fmt_half <- function(x) {
    if (is.null(digits)) return(x / 2)
    sub("\\\\s+$", "", formatC(x / 2, digits = digits, format = "g"))
  }

  # Remove the token body from `txt` and consume one adjacent comma if present.
  # span = c(start, end) (1-based inclusive indices) of the token body **without** the leading delimiter.
  splice_out_token <- function(txt, span) {
    s <- span[1]; e <- span[2]; n <- nchar(txt)
    ch_at <- function(i) if (i < 1 || i > n) "" else substr(txt, i, i)
    left  <- ch_at(s - 1)
    right <- ch_at(e + 1)

    if (left == ",") {
      # Prefer consuming the comma BEFORE the token: [... , <token> ...] -> [...]
      paste0(substr(txt, 1, s - 2), substr(txt, e + 1, n))
    } else if (right == ",") {
      # Otherwise consume the comma AFTER the token: [... <token> , ...] -> [...]
      paste0(substr(txt, 1, s - 1), substr(txt, e + 2, n))
    } else {
      # No adjacent comma: just excise the token body itself.
      paste0(substr(txt, 1, s - 1), substr(txt, e + 1, n))
    }
  }

  # Find the outgroup token in ONE of two strict forms and return span + L:
  #   a) "(label:number"  -> body starts right after "("
  #   b) ",label:number"  -> body starts right after comma (NO SPACES allowed)
  # We explicitly DO NOT include the leading delimiter inside the span.
  find_outgroup_token <- function(txt, label) {
    lab_re <- re_escape(label)

    # --- Case a) "(label:number"
    pat_a <- paste0("\\(", lab_re, ":[0-9Ee+\\.-]+")
    ma <- regexpr(pat_a, txt, perl = TRUE)
    if (ma[1] != -1L) {
      start <- as.integer(ma[1])
      len   <- attr(ma, "match.length")[1]
      end   <- start + len - 1L
      span  <- c(start + 1L, end)  # skip '('
      num   <- suppressWarnings(as.numeric(sub(".*:([0-9Ee+\\.-]+)\\s*$", "\\1",
                                               substr(txt, start, end), perl = TRUE)))
      if (is.finite(num) && num >= 0) return(list(span = span, L = num))
    }

    # --- Case b) ",label:number"  (STRICT: no spaces after comma)
    pat_b <- paste0(",", lab_re, ":[0-9Ee+\\.-]+")
    mb <- regexpr(pat_b, txt, perl = TRUE)
    if (mb[1] != -1L) {
      start <- as.integer(mb[1])
      len   <- attr(mb, "match.length")[1]
      end   <- start + len - 1L
      span  <- c(start + 1L, end)  # skip ',' exactly one char; no spaces are allowed
      num   <- suppressWarnings(as.numeric(sub(".*:([0-9Ee+\\.-]+)\\s*$", "\\1",
                                               substr(txt, start, end), perl = TRUE)))
      if (is.finite(num) && num >= 0) return(list(span = span, L = num))
    }

    NULL
  }

  ## --------------------------- 1) Read input safely ---------------------------
  if (!file.exists(iqtree_treefile)) stop("Input file does not exist.")
  tree_txt <- read_one_line(iqtree_treefile)
  # Normalize: ensure a trailing semicolon for predictable downstream ops.
  if (!grepl(";\\s*$", tree_txt)) tree_txt <- paste0(gsub("\\s+$", "", tree_txt), ";")

  ## --------------------------- 2) Locate token & L ----------------------------
  hit <- find_outgroup_token(tree_txt, outgroup_label)
  if (is.null(hit)) {
    stop(
      "Outgroup label not found in strict forms. Expecting either ",
      "'(", outgroup_label, ":<len>' or ',", outgroup_label, ":<len>' (no spaces after comma)."
    )
  }
  L <- hit$L

  ## --------------------------- 3) Remove original token -----------------------
  ingroup_txt <- splice_out_token(tree_txt, hit$span)
  ingroup_txt <- sub(";\\s*$", "", ingroup_txt)  # drop terminal semicolon

  ## --------------------------- 4) Build new rooted tree -----------------------
  half <- fmt_half(L)
  newick <- paste0("(", ingroup_txt, ",", outgroup_label, ":", half, "):", half, ";")

  ## --------------------------- 5) Write output --------------------------------
  outfp <- sub("\\.treefile$", ".bi.treefile", iqtree_treefile)
  if (identical(outfp, iqtree_treefile)) outfp <- paste0(iqtree_treefile, ".bi.treefile")
  writeLines(newick, outfp)
  invisible(outfp)
}
