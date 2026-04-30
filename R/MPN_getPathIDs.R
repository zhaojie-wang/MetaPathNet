#' Retrieve and Classify Organism-Specific KEGG Pathways
#'
#' Retrieves KEGG pathway IDs for one or more organisms and assigns each to a
#' functional category and type based on KEGG BRITE hierarchy.
#'
#' @param organism_code Character vector of KEGG organism codes.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{Path_id}{Full KEGG pathway ID.}
#'   \item{Path_description}{Full KEGG pathway name.}
#'   \item{Path_category}{Top-level and subcategory from KEGG BRITE hierarchy.}
#'   \item{Path_type}{\code{"metabolic"} or \code{"signaling"}, based on top-level category.}
#'   \item{Species}{KEGG organism code.}
#' }
#'
#' @details
#' Pathway classification uses the KEGG BRITE hierarchy \code{br08901}.
#' Pathways under the top-level category \code{"Metabolism"} are typed as
#' \code{"metabolic"}; all others are typed as \code{"signaling"}. Global
#' and overview maps are excluded from the output.
#'
#' @examples
#' ## 1) Retrieve all human pathways and filter for tryptophan metabolism
#' path_hsa <- MPN_getPathIDs("hsa")
#' trp_hsa <- subset(path_hsa, grepl("Tryptophan metabolism", Path_description))
#' head(trp_hsa)
#'
#' ## 2) Retrieve pathways for multiple organisms
#' org_vec <- c("hsa", "eco", "bsu")
#' path_multi <- MPN_getPathIDs(org_vec)
#'
#' ## 3) Show pathways shared across the selected organisms
#' ##    (appear once per species, so shared ones have count == length(org_vec))
#' shared_counts <- table(path_multi$Path_description)
#' head(shared_counts[shared_counts == length(org_vec)])
#'
#' ## 4) Filter tryptophan metabolism pathways and show which species contain them
#' trp_multi <- subset(path_multi, grepl("Tryptophan metabolism", Path_description))
#' head(trp_multi)
#' table(trp_multi$Species)
#'
#' @seealso
#' \code{\link{MPN_keggNetwork}} for building networks from these pathway sets.
#' @importFrom KEGGREST keggGet keggLink
#' @export
MPN_getPathIDs <- function(organism_code) {
  organism_code <- unique(as.character(organism_code))

  #==============================================================
  # STEP 1 — Parse KEGG BRITE (br08901) into a pathway catalogue
  #==============================================================
  # br08901 holds the global pathway hierarchy (Metabolism vs others).
  brite_raw <- tryCatch(
    KEGGREST::keggGet("br:br08901")[[1]],
    error = function(e) stop("Failed to retrieve KEGG BRITE hierarchy (br08901).")
  )

  # Normalise to a clean character vector of lines
  lines <- if (length(brite_raw) == 1L) {
    strsplit(brite_raw, "\n", fixed = TRUE)[[1]]
  } else {
    brite_raw
  }

  lines <- sub("\r$", "", lines)
  trim  <- sub("^\\s+", "", lines)

  # Track current A (top level) and B (sub-category) for each line
  current_A <- current_B <- character(length(trim))
  lastA <- lastB <- NA_character_

  for (i in seq_along(trim)) {
    if (grepl("^A", trim[i])) {
      lastA <- sub("^A", "", trim[i])
    } else if (grepl("^B", trim[i])) {
      lastB <- sub("^B\\s+", "", trim[i])
    }
    current_A[i] <- lastA
    current_B[i] <- lastB
  }

  # Keep only C-level lines
  is_C    <- grepl("^C\\s+[0-9]{5}\\b", trim)
  C_lines <- trim[is_C]

  # Extract ID, description, and coarse type (metabolic / signaling)
  Path_id   <- sub("^C\\s+([0-9]{5}).*$", "\\1", C_lines)
  Path_desc <- sub("^C\\s+[0-9]{5}\\s+", "", C_lines)
  Path_type <- ifelse(current_A[is_C] == "Metabolism", "metabolic", "signaling")

  # Merge A and B levels into a single category string
  Path_cat <- paste(current_A[is_C], current_B[is_C], sep = "; ")

  path_info_df <- data.frame(
    Path_id          = as.character(Path_id),
    Path_description = Path_desc,
    Path_category    = Path_cat,
    Path_type        = Path_type,
    stringsAsFactors = FALSE
  )

  # Drop global / overview maps (keep only “real” functional maps)
  path_info_df <- path_info_df[
    !grepl("Global and overview maps", path_info_df$Path_category),
  ]

  # Add suffix (numeric pathway code) for joining with organism-specific IDs
  general_info        <- path_info_df
  general_info$Suffix <- general_info$Path_id

  #==============================================================
  # STEP 2 — Attach organism-specific pathway IDs (keggLink)
  #==============================================================
  all_species_paths <- lapply(organism_code, function(org) {

    # Pull all pathway links for this organism
    sp_link <- tryCatch(
      KEGGREST::keggLink("pathway", org),
      error = function(e) NULL
    )

    if (is.null(sp_link) || length(sp_link) == 0L) {
      message("[", org, "] Retrieved 0 pathway IDs.")
      return(NULL)
    }

    sp_paths <- unique(gsub("path:", "", sp_link))
    message("[", org, "] Retrieved ", length(sp_paths), " pathway IDs.")

    suffixes <- sub("^[a-zA-Z]+", "", sp_paths)
    df <- data.frame(
      Species_Path_id = sp_paths,
      Suffix          = suffixes,
      Species         = org,
      stringsAsFactors = FALSE
    )

    # Join with the global pathway catalogue on numeric suffix
    merged <- merge(df, general_info, by = "Suffix", all.x = TRUE)
    merged <- merged[!is.na(merged$Path_description), ]

    if (nrow(merged) == 0L) {
      return(NULL)
    }

    merged <- merged[, c(
      "Species_Path_id", "Path_description", "Path_category", "Path_type", "Species"
    )]
    colnames(merged)[1] <- "Path_id"

    merged
  })

  #==============================================================
  # STEP 3 — Combine species tables and tidy output
  #==============================================================
  ## Summary message across organisms
  ok_idx   <- !vapply(all_species_paths, is.null, logical(1))
  ok_orgs  <- organism_code[ok_idx]
  bad_orgs <- organism_code[!ok_idx]

  message(
    "Pathway retrieval summary: ",
    length(ok_orgs), "/", length(organism_code), " organism(s) succeeded",
    if (length(bad_orgs) > 0) paste0(" (failed: ", paste(bad_orgs, collapse = ", "), ")") else ""
  )

  all_species_paths <- Filter(Negate(is.null), all_species_paths)

  if (length(all_species_paths) == 0L) {
    stop("No valid KEGG pathways could be retrieved for the supplied organism code(s).")
  }

  final_df <- do.call(rbind, all_species_paths)
  rownames(final_df) <- NULL
  final_df[] <- lapply(final_df, as.character)

  return(final_df)
}

