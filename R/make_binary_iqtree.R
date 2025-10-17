make_binary_iqtree <- function(iqtree_treefile, outgroup_label, digits = NULL) {
  # ---------- tiny helpers ----------
  read_one_line <- function(p) paste(readLines(p, warn = FALSE), collapse = "")
  fmt_half <- function(x) {
    if (is.null(digits)) return(x/2)
    sub("\\s+$", "", formatC(x/2, digits = digits, format = "g"))
  }
  # escape regex metacharacters in label (safer literal match)
  re_escape <- function(x) gsub("([][{}()+*^$|?\\.^-])", "\\\\\\\\1", x, perl = TRUE)

  # ---------- read & normalize ----------
  if (!file.exists(iqtree_treefile)) stop("Input file does not exist: ", iqtree_treefile)
  txt <- read_one_line(iqtree_treefile)
  if (!grepl(";\\s*$", txt)) txt <- paste0(gsub("\\s+$", "", txt), ";")

  # ---------- strict patterns (no spaces after comma) ----------
  lab <- re_escape(outgroup_label)
  num <- "([0-9]+(?:\\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)"

  # A: "(label:length,"     -> remove whole token and following comma, keep "("
  patA <- paste0("\\(", lab, ":", num, ",")
  # B: ",label:length);"     -> remove preceding comma + token, keep ")"
  patB <- paste0(",",  lab, ":", num, "\\)\\s*;\\s*$")
  # C: ",label:length,"      -> remove preceding comma + token + following comma; keep a single ","
  patC <- paste0(",",  lab, ":", num, ",")

  L <- NA_real_
  ingroup <- NULL

  mA <- regexec(patA, txt, perl = TRUE)
  if (mA[[1]][1] != -1L) {
    cap <- regmatches(txt, mA)[[1]]
    L <- as.numeric(cap[2])
    ingroup <- sub(patA, "(", txt, perl = TRUE)
  } else {
    mB <- regexec(patB, txt, perl = TRUE)
    if (mB[[1]][1] != -1L) {
      cap <- regmatches(txt, mB)[[1]]
      L <- as.numeric(cap[2])
      ingroup <- sub(patB, ")", txt, perl = TRUE)
    } else {
      mC <- regexec(patC, txt, perl = TRUE)
      if (mC[[1]][1] != -1L) {
        cap <- regmatches(txt, mC)[[1]]
        L <- as.numeric(cap[2])
        # replace ",label:length," with a single comma to bridge neighbors
        ingroup <- sub(patC, ",", txt, perl = TRUE)
      }
    }
  }

  if (!is.finite(L)) {
    stop(paste0(
      "Outgroup not found in strict forms.\n",
      "Expected one of:\n",
      "  (", outgroup_label, ":<len>,\n",
      "  ,", outgroup_label, ":<len>);\n",
      "  ,", outgroup_label, ":<len>,\n",
      "with NO spaces after commas."
    ))
  }

  # drop trailing semicolon from ingroup
  ingroup <- sub(";\\s*$", "", ingroup)

  # ---------- build final binary-rooted tree ----------
  half <- fmt_half(L)
  newick <- paste0("((", ingroup, "):", half, ",", outgroup_label, ":", half, ");")

  # ---------- write output ----------
  outfp <- sub("\\.treefile$", ".bi.treefile", iqtree_treefile)
  if (identical(outfp, iqtree_treefile)) outfp <- paste0(iqtree_treefile, ".bi.treefile")
  writeLines(newick, outfp)
  invisible(outfp)
}
