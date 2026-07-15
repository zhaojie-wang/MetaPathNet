#' Suggest Candidate Prokaryotic Organisms or KOs from KEGG Identifiers
#'
#' Suggest candidate prokaryotic organisms from KEGG KO IDs or candidate
#' KEGG Orthology (KO) annotations from KEGG compound IDs.
#'
#' The organism-suggestion workflow focuses on prokaryotic organisms to
#' support host-microbiome network construction. \emph{Homo sapiens}
#' (\code{hsa}) is retained when present as the host reference.
#'
#' @param query Character vector of KEGG identifiers. For
#'   \code{entity = "ko"}, KO IDs should follow the form
#'   \code{"K#####"}. For \code{entity = "compound"}, compound IDs
#'   should follow the form \code{"C#####"} or \code{"cpd:C#####"}.
#' @param entity Character string specifying the input type:
#'   \code{"ko"} or \code{"compound"}. Only one entity type is
#'   accepted per call.
#'
#' @return A \code{data.frame} whose structure depends on
#'   \code{entity}:
#' \describe{
#'   \item{For \code{entity = "ko"}:}{
#'     Columns \code{organism_name}, \code{organism_code},
#'     \code{overlap_count}, and \code{overlap_hits}. Candidate
#'     prokaryotic organisms and \code{hsa}, when present, are ranked
#'     by decreasing \code{overlap_count}.
#'   }
#'   \item{For \code{entity = "compound"}:}{
#'     Columns \code{compound_id}, \code{compound_name},
#'     \code{KO_ID}, and \code{KO_name}. Duplicate rows are removed.
#'   }
#' }
#' An empty data frame with the appropriate columns is returned when
#' no candidates are found or when the required KEGG genome
#' classification cannot be retrieved.
#'
#' @details
#' For \code{entity = "ko"}, each KO is queried with
#' \code{KEGGREST::keggLink()} to retrieve genes carrying that
#' orthology assignment. Organism codes are extracted from the linked
#' KEGG gene IDs and matched to the live KEGG genome table.
#'
#' Broad organism classification is reconstructed from the KEGG
#' eukaryotic and prokaryotic genome groups. Candidate organisms are
#' restricted to prokaryotes, while \emph{Homo sapiens}
#' (\code{hsa}) is retained as a host reference. Organisms are ranked
#' according to the number of input KOs represented in each genome.
#'
#' For \code{entity = "compound"}, each compound is queried with
#' \code{KEGGREST::keggGet()}. Linked KEGG reaction IDs are extracted,
#' and the \code{ORTHOLOGY} field of each reaction is used to retrieve
#' associated candidate KOs.
#'
#' Input IDs are validated before querying KEGG. If none of the supplied
#' IDs has a valid format, the function stops with an error.
#'
#' @examples
#' \donttest{
#' ## 1) Retrieve candidate KOs associated with selected compounds
#' cpd_ids <- MPN_keggFinder(
#'   KEGG_database = "compound",
#'   searchBy      = "name",
#'   query         = c("tryptophan", "indole", "serotonin")
#' )
#'
#' ko_results <- MPN_suggestEntities(
#'   query  = cpd_ids$KEGG_ID,
#'   entity = "compound"
#' )
#'
#' head(ko_results)
#'
#' ## 2) Retrieve candidate prokaryotic organisms from selected KOs
#' org_results <- MPN_suggestEntities(
#'   query  = c("K18277", "K07811", "K00108", "K14156"),
#'   entity = "ko"
#' )
#'
#' head(org_results)
#' }
#'
#' @seealso \code{\link{MPN_keggFinder}} for converting compound names,
#'   organism names, or enzyme EC numbers to KEGG identifiers.
#'
#' @importFrom KEGGREST keggGet keggLink keggList
#' @export
MPN_suggestEntities <- function(
    query,
    entity = c("ko", "compound")
) {

  #==============================================================
  # Step 1 - Validate 'entity' argument
  #==============================================================
  if (missing(entity) || length(entity) != 1L || !is.character(entity)) {
    stop(
      "'entity' must be a single character string: either \"ko\" or \"compound\".\n",
      "Only one entity type is accepted per call - merging is not allowed."
    )
  }

  entity <- tolower(trimws(entity))

  if (!entity %in% c("ko", "compound")) {
    stop(
      "Unrecognised entity type: \"", entity, "\".\n",
      "Accepted values are \"ko\" (to find candidate organisms) ",
      "or \"compound\" (to find candidate KOs)."
    )
  }

  #==============================================================
  # Step 2 - Validate and normalise 'query' input
  #==============================================================
  if (missing(query) || length(query) == 0L) {
    stop("'query' must be a non-empty character vector of KEGG IDs.")
  }

  query <- unique(trimws(as.character(query)))
  query <- query[!is.na(query) & nzchar(query)]

  if (length(query) == 0L) {
    stop("All provided query IDs are missing or empty after trimming.")
  }

  ## Entity-specific format check
  if (entity == "ko") {
    valid_mask <- grepl("^K\\d{5}$", query)
    fmt_hint   <- "KO IDs must match the pattern 'K#####' (e.g. K00169)."
  } else {
    valid_mask <- grepl("^(cpd:)?[Cc][0-9]{5}$", query)
    fmt_hint   <- "Compound IDs must match 'C#####' or 'cpd:C#####' (e.g. C00114)."
  }

  if (!any(valid_mask)) {
    stop(
      "None of the supplied IDs have a valid format for entity=\"", entity, "\".\n",
      fmt_hint, "\n",
      "First problematic IDs: ",
      paste(head(query[!valid_mask], 5), collapse = ", ")
    )
  }

  unmatched_format <- query[!valid_mask]
  query <- query[valid_mask]

  if (length(unmatched_format) > 0L) {
    message(
      "Note: ", length(unmatched_format),
      " ID(s) removed due to invalid format:\n  ",
      paste(unmatched_format, collapse = ", "), "\n",
      fmt_hint
    )
  }

  # ========================  entity == "ko"  ======================
  if (entity == "ko") {

    empty_result <- data.frame(
      organism_name = character(0),
      organism_code = character(0),
      overlap_count = integer(0),
      overlap_hits  = character(0),
      stringsAsFactors = FALSE
    )

    #==============================================================
    # Step 3 - Retrieve the live KEGG genome tables
    #==============================================================
    genome <- tryCatch(
      KEGGREST::keggList("genome"),
      error = function(e) NULL
    )

    prokaryotes <- tryCatch(
      KEGGREST::keggList("genome", "prokaryotes"),
      error = function(e) NULL
    )

    eukaryotes <- tryCatch(
      KEGGREST::keggList("genome", "eukaryotes"),
      error = function(e) NULL
    )

    if (
      is.null(genome) ||
      length(genome) == 0L ||
      is.null(names(genome)) ||
      is.null(prokaryotes) ||
      length(prokaryotes) == 0L ||
      is.null(names(prokaryotes)) ||
      is.null(eukaryotes) ||
      length(eukaryotes) == 0L ||
      is.null(names(eukaryotes))
    ) {
      warning(
        "KEGG genome classification could not be retrieved; ",
        "returning an empty organism candidate table."
      )
      return(empty_result)
    }

    #==============================================================
    # Step 4 - Build internal organism table
    #==============================================================
    genome_value <- trimws(as.character(genome))

    valid_genome <- !is.na(names(genome)) &
      nzchar(names(genome)) &
      !is.na(genome_value) &
      nzchar(genome_value) &
      grepl(";", genome_value, fixed = TRUE)

    if (!any(valid_genome)) {
      warning(
        "KEGG genome table returned no parseable entries; ",
        "returning an empty organism candidate table."
      )
      return(empty_result)
    }

    genome <- genome[valid_genome]
    genome_value <- genome_value[valid_genome]

    ## Reconstruct broad genome-group classification from KEGG genome groups
    organism_df <- data.frame(
      taxon_id = names(genome),
      organism = trimws(sub(";.*$", "", genome_value)),
      species = trimws(sub("^[^;]+;\\s*", "", genome_value)),
      phylogeny = NA_character_,
      stringsAsFactors = FALSE
    )

    organism_df <- organism_df[
      nzchar(organism_df$taxon_id) &
        nzchar(organism_df$organism) &
        nzchar(organism_df$species),
      ,
      drop = FALSE
    ]

    if (nrow(organism_df) == 0L) {
      warning(
        "KEGG genome table returned no valid organism entries; ",
        "returning an empty organism candidate table."
      )
      return(empty_result)
    }

    prokaryote_value <- trimws(as.character(prokaryotes))

    valid_prokaryote <- !is.na(names(prokaryotes)) &
      nzchar(names(prokaryotes)) &
      !is.na(prokaryote_value) &
      nzchar(prokaryote_value) &
      grepl(";", prokaryote_value, fixed = TRUE)

    eukaryote_value <- trimws(as.character(eukaryotes))

    valid_eukaryote <- !is.na(names(eukaryotes)) &
      nzchar(names(eukaryotes)) &
      !is.na(eukaryote_value) &
      nzchar(eukaryote_value) &
      grepl(";", eukaryote_value, fixed = TRUE)

    if (!any(valid_prokaryote) || !any(valid_eukaryote)) {
      warning(
        "KEGG genome groups returned no parseable entries; ",
        "returning an empty organism candidate table."
      )
      return(empty_result)
    }

    prokaryote_codes <- unique(trimws(
      sub(";.*$", "", prokaryote_value[valid_prokaryote])
    ))

    eukaryote_codes <- unique(trimws(
      sub(";.*$", "", eukaryote_value[valid_eukaryote])
    ))

    organism_df$phylogeny[
      organism_df$organism %in% eukaryote_codes
    ] <- "Eukaryotes;"

    organism_df$phylogeny[
      organism_df$organism %in% prokaryote_codes
    ] <- "Prokaryotes;"

    rownames(organism_df) <- NULL

    ## Remove parenthesised common names from species labels
    has_parentheses <- grepl("\\(.*\\)", organism_df$species)

    organism_df$species[has_parentheses] <- trimws(
      sub(
        "\\s*\\(.*\\)$",
        "",
        organism_df$species[has_parentheses]
      )
    )

    #==============================================================
    # Step 5 - Query KEGG for organism codes linked to each KO
    #==============================================================
    ko_org_list <- vector("list", length(query))
    names(ko_org_list) <- query
    failed_kos <- character(0)

    for (i in seq_along(query)) {
      ko <- query[i]

      ## Respect KEGG REST API rate limits for larger input sets
      if (length(query) > 10L && i > 1L && (i - 1L) %% 10L == 0L) {
        Sys.sleep(1)
      }

      ko_data <- tryCatch(
        KEGGREST::keggLink("genes", ko),
        error = function(e) NULL
      )

      if (is.null(ko_data)) {
        failed_kos <- c(failed_kos, ko)
        ko_org_list[[ko]] <- character(0)
        next
      }

      if (length(ko_data) == 0L) {
        ko_org_list[[ko]] <- character(0)
        next
      }

      ## Extract organism codes from gene IDs
      gene_ids <- as.character(ko_data)
      org_codes <- sub(":.*$", "", gene_ids)

      ko_org_list[[ko]] <- unique(
        org_codes[!is.na(org_codes) & nzchar(org_codes)]
      )
    }

    ## Report retrieval failures
    if (length(failed_kos) > 0L) {
      message(
        "KEGG query failed for ", length(failed_kos), " KO(s): ",
        paste(failed_kos, collapse = ", "),
        "\nThese were skipped."
      )
    }

    #==============================================================
    # Step 6 - Identify KOs with zero organism hits
    #==============================================================
    empty_kos <- names(ko_org_list)[
      vapply(ko_org_list, length, integer(1)) == 0L
    ]

    if (length(empty_kos) > 0L) {
      message(
        "Note: ", length(empty_kos),
        " KO(s) returned no organism links:\n  ",
        paste(empty_kos, collapse = ", ")
      )
    }

    #==============================================================
    # Step 7 - Collect all unique candidate organism codes
    #==============================================================
    all_org_codes <- unique(unlist(
      ko_org_list,
      use.names = FALSE
    ))

    if (length(all_org_codes) == 0L) {
      warning(
        "No organism candidates could be retrieved for any ",
        "of the provided KO IDs."
      )
      return(empty_result)
    }

    #==============================================================
    # Step 8 - Build overlap summary
    #==============================================================
    out_df <- do.call(
      rbind,
      lapply(all_org_codes, function(org_code) {
        matched_kos <- names(ko_org_list)[vapply(
          ko_org_list,
          function(x) org_code %in% x,
          logical(1)
        )]

        data.frame(
          organism_code = org_code,
          overlap_count = length(matched_kos),
          overlap_hits = paste(matched_kos, collapse = ";"),
          stringsAsFactors = FALSE
        )
      })
    )

    #==============================================================
    # Step 9 - Attach species names and retain prokaryotes only
    #==============================================================
    match_idx <- match(
      out_df$organism_code,
      organism_df$organism
    )

    out_df$organism_name <- ifelse(
      !is.na(match_idx) &
        !is.na(organism_df$species[match_idx]) &
        nzchar(organism_df$species[match_idx]),
      organism_df$species[match_idx],
      out_df$organism_code
    )

    out_df$phylogeny <- organism_df$phylogeny[match_idx]

    ## Filter to prokaryotes and Homo sapiens
    out_df <- out_df[
      (
        !is.na(out_df$phylogeny) &
          grepl("^Prokaryotes;", out_df$phylogeny)
      ) |
        out_df$organism_code == "hsa",
      ,
      drop = FALSE
    ]

    if (nrow(out_df) == 0L) {
      message(
        "No prokaryotic organisms or Homo sapiens candidates remained ",
        "after genome-group filtering."
      )
      return(empty_result)
    }

    #==============================================================
    # Step 10 - Final formatting and return
    #==============================================================
    out_df <- out_df[
      ,
      c(
        "organism_name",
        "organism_code",
        "overlap_count",
        "overlap_hits"
      ),
      drop = FALSE
    ]

    out_df <- out_df[
      order(-out_df$overlap_count, out_df$organism_code),
      ,
      drop = FALSE
    ]

    rownames(out_df) <- NULL

    return(out_df)
  }

  # =====================  entity == "compound"  ===================

  #==============================================================
  # Step 3 - Standardise compound IDs to "cpd:CXXXXX"
  #==============================================================
  query <- toupper(gsub("^cpd:", "", query))
  query <- paste0("cpd:", query)

  #==============================================================
  # Step 4 - Loop over compounds -> reactions -> KOs
  #==============================================================
  out_rows <- list()
  row_index <- 1L
  failed_cpds <- character(0)
  no_reaction <- character(0)

  for (i in seq_along(query)) {
    cpd <- query[i]

    ## Respect KEGG REST API rate limits for larger input sets
    if (length(query) > 10L && i > 1L && (i - 1L) %% 10L == 0L) {
      Sys.sleep(1)
    }

    ## Retrieve compound entry
    compound_info <- tryCatch(
      KEGGREST::keggGet(cpd),
      error = function(e) NULL
    )

    if (
      is.null(compound_info) ||
      length(compound_info) == 0L ||
      is.null(compound_info[[1]])
    ) {
      failed_cpds <- c(failed_cpds, cpd)
      next
    }

    ## Extract compound name
    compound_name <- cpd

    if (
      !is.null(compound_info[[1]]$NAME) &&
      length(compound_info[[1]]$NAME) > 0L
    ) {
      compound_name <- sub(
        ";.*$",
        "",
        compound_info[[1]]$NAME[1]
      )
    }

    ## Extract linked reaction IDs
    if (
      is.null(compound_info[[1]]$REACTION) ||
      length(compound_info[[1]]$REACTION) == 0L
    ) {
      no_reaction <- c(no_reaction, cpd)
      next
    }

    reaction_ids <- unique(trimws(unlist(
      strsplit(compound_info[[1]]$REACTION, " ")
    )))

    reaction_ids <- reaction_ids[
      !is.na(reaction_ids) & nzchar(reaction_ids)
    ]

    if (length(reaction_ids) == 0L) {
      no_reaction <- c(no_reaction, cpd)
      next
    }

    ## Collect orthology entries for each linked reaction
    for (j in seq_along(reaction_ids)) {
      rid <- reaction_ids[j]

      if (
        length(reaction_ids) > 10L &&
        j > 1L &&
        (j - 1L) %% 10L == 0L
      ) {
        Sys.sleep(1)
      }

      reaction_info <- tryCatch(
        KEGGREST::keggGet(rid),
        error = function(e) NULL
      )

      if (
        is.null(reaction_info) ||
        length(reaction_info) == 0L ||
        is.null(reaction_info[[1]])
      ) {
        next
      }

      if (
        is.null(reaction_info[[1]]$ORTHOLOGY) ||
        length(reaction_info[[1]]$ORTHOLOGY) == 0L
      ) {
        next
      }

      ko_ids <- names(reaction_info[[1]]$ORTHOLOGY)
      ko_names <- as.character(reaction_info[[1]]$ORTHOLOGY)

      if (length(ko_ids) == 0L) {
        next
      }

      out_rows[[row_index]] <- data.frame(
        compound_id = cpd,
        compound_name = compound_name,
        KO_ID = ko_ids,
        KO_name = ko_names,
        stringsAsFactors = FALSE
      )

      row_index <- row_index + 1L
    }
  }

  #==============================================================
  # Step 5 - Report retrieval failures and unmatched compounds
  #==============================================================
  if (length(failed_cpds) > 0L) {
    message(
      "KEGG query failed for ", length(failed_cpds),
      " compound(s): ",
      paste(failed_cpds, collapse = ", "),
      "\nThese were skipped."
    )
  }

  if (length(no_reaction) > 0L) {
    message(
      "Note: ", length(no_reaction),
      " compound(s) had no linked reactions:\n  ",
      paste(no_reaction, collapse = ", ")
    )
  }

  if (length(out_rows) == 0L) {
    out_df <- data.frame(
      compound_id = character(0),
      compound_name = character(0),
      KO_ID = character(0),
      KO_name = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    out_df <- do.call(rbind, out_rows)
  }

  ## Identify compounds with no associated candidate KOs
  matched_cpds <- unique(out_df$compound_id)
  unmatched_cpds <- setdiff(query, matched_cpds)

  if (length(unmatched_cpds) > 0L) {
    message(
      "Summary: ", length(unmatched_cpds), " of ", length(query),
      " compound(s) had no associated candidate KOs:\n  ",
      paste(unmatched_cpds, collapse = ", ")
    )
  }

  #==============================================================
  # Step 6 - Deduplicate and return
  #==============================================================
  if (nrow(out_df) == 0L) {
    warning(
      "No candidate KOs could be retrieved for any ",
      "of the provided compound IDs."
    )
    return(out_df)
  }

  out_df <- unique(out_df)
  rownames(out_df) <- NULL

  out_df
}
