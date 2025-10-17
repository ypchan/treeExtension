make_binary_iqtree <- function(tre_file, outgroup_label) {
    # 1. Read the tree file into a single long string
    tre_string <- readLines(tre_file, warn = FALSE)
    tre_string <- paste(tre_string, collapse = "")  # Combine into one-line string
    
    # Check if outgroup_label exists in the tree string
    if (!grepl(outgroup_label, tre_string, fixed = TRUE)) {
        stop(paste(outgroup_label, "not in tre_file"))
    }
    
    # 2. Extract pattern and convert L to numeric
    # Escape special characters in outgroup label for regex compatibility
    escaped_label <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outgroup_label)
    pattern <- paste0(escaped_label, ":0\\.[0-9]+")
    match_result <- regexpr(pattern, tre_string, perl = TRUE)
    
    if (match_result == -1) {
        stop(paste("No matching pattern found:", pattern))
    }
    
    full_match <- regmatches(tre_string, match_result)
    L_str <- sub(paste0(outgroup_label, ":"), "", full_match)  # Extract numeric portion as string
    L <- as.numeric(L_str)  # Convert L to numeric value
    
    if (is.na(L)) {
        stop("Extracted L could not be converted to a numeric value")
    }
    
    # 3. Attempt to remove specified strings in sequence
    target1 <- paste0(",", outgroup_label, ":", L_str)  # Use original string for exact matching
    target2 <- paste0(outgroup_label, ":", L_str, ",")
    
    if (grepl(target1, tre_string, fixed = TRUE)) {
        tre_string <- gsub(target1, "", tre_string, fixed = TRUE)
    } else if (grepl(target2, tre_string, fixed = TRUE)) {
        tre_string <- gsub(target2, "", tre_string, fixed = TRUE)
    } else {
        stop(paste("No removable target strings found:", target1, "or", target2))
    }
    
    # 4. Remove semicolon and record as ingroup_string
    ingroup_string <- gsub(";", "", tre_string, fixed = TRUE)
    
    # 5. Combine to form the binary tree string with numeric L/2
    L_half <- L / 2  # Calculate half of L as numeric
    binary_tre_string <- paste0("(", ingroup_string, ":", L_half, ",", outgroup_label, ":", L_half, ");")
    
    # 6. Create output file path
    input_basename <- basename(tre_file)
    output_filename <- paste0(tools::file_path_sans_ext(input_basename), ".bi.treefile")
    output_path <- file.path(dirname(tre_file), output_filename)
    
    # 7. Write result to output file
    writeLines(binary_tre_string, con = output_path)
    
    message(paste("Output written to:", output_path))
}
