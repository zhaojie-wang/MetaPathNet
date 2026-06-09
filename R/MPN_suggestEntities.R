#' Suggest Candidate Organisms or KOs from KEGG Identifiers
#'
#' Suggest candidate organisms from KEGG KO IDs or candidate KEGG
#' Orthology (KO) annotations from KEGG compound IDs.
#'
#' This function supports semi-automated candidate selection during
#' network construction.
#'
#' @param query Character vector of KEGG identifiers. For
#'   \code{entity = "ko"}, KO IDs should be in the form
#'   \code{"K#####"}. For \code{entity = "compound"}, compound IDs
#'   should be in the form \code{"C#####"} or \code{"cpd:C#####"}.
#' @param entity Character string specifying the input type:
#'   \code{"ko"} or \code{"compound"}. Only one entity type is
#'   accepted per call.
#'
#' @return A \code{data.frame} whose structure depends on the entity type:
#' \describe{
#'   \item{For \code{entity = "ko"}:}{Columns \code{organism_name},
#'     \code{organism_code}, \code{overlap_count}, and
#'     \code{overlap_hits}. Rows are ranked by decreasing
#'     \code{overlap_count}.}
#'   \item{For \code{entity = "compound"}:}{Columns
#'     \code{compound_id}, \code{compound_name}, \code{KO_ID}, and
#'     \code{KO_name}. Rows are deduplicated.}
#' }
#' An empty data frame with the appropriate columns is returned when
#' no candidates are found.
#'
#' @details
#' For \code{entity = "ko"}, each KO is queried with
#' \code{KEGGREST::keggLink()} to retrieve genes carrying that
#' orthology assignment. Organism codes are extracted from the linked
#' KEGG gene IDs, matched to the live KEGG organism table, filtered to
#' prokaryotes and \emph{Homo sapiens} (\code{hsa}), and ranked by the
#' number of query KOs found in each organism.
#'
#' For \code{entity = "compound"}, each compound is queried with
#' \code{KEGGREST::keggGet()}, linked KEGG reaction IDs are extracted,
#' and each reaction is searched for KO annotations through its
#' \code{$ORTHOLOGY} field.
#'
#' Input IDs are validated before querying KEGG. If none of the input
#' IDs have a valid format, the function stops with an error. Invalid
#' IDs are otherwise removed with a message. KEGG retrieval failures,
#' missing reaction links, and missing orthology annotations are
#' handled internally and reported without interrupting execution.
#'
#' @examples
#' \donttest{
#' ## 1) Compound to KO: retrieve candidate KOs from mapped compounds
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
#' ## 2) KO to organism: retrieve candidate prokaryotic organisms
#' org_results <- MPN_suggestEntities(
#'   query  = ko_results$KO_ID,
#'   entity = "ko"
#' )
#'
#' head(org_results)
#' }
#'
#' @seealso \code{\link{MPN_keggFinder}} for converting metabolite names,
#'   organism names, or enzyme EC numbers to KEGG identifiers.
#'
#' @importFrom KEGGREST keggGet keggLink
#' @export
MPN_suggestEntities <- function(query, entity = c("ko", "compound")) {

  #==============================================================
  # Step 1 - Validate 'entity' argument
  #==============================================================
  if (missing(entity) || length(entity) != 1 || !is.character(entity)) {
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
  if (missing(query) || length(query) == 0) {
    stop("'query' must be a non-empty character vector of KEGG IDs.")
  }

  query <- unique(trimws(as.character(query)))
  query <- query[query != ""]

  if (length(query) == 0) {
    stop("All provided query IDs are empty after trimming.")
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
      "First problematic IDs: ", paste(head(query[!valid_mask], 5), collapse = ", ")
    )
  }

  unmatched_format <- query[!valid_mask]
  query <- query[valid_mask]

  if (length(unmatched_format) > 0) {
    message(
      "Note: ", length(unmatched_format), " ID(s) removed due to invalid format:\n  ",
      paste(unmatched_format, collapse = ", "), "\n",
      fmt_hint
    )
  }

  # ========================  entity == "ko"  ======================
  if (entity == "ko") {

    #==============================================================
    # Step 3 - Retrieve the live KEGG organism table
    #==============================================================
    response <- tryCatch(
      rawToChar(
        curl::curl_fetch_memory(
          "https://rest.kegg.jp/list/organism",
          handle = curl::new_handle()
        )$content
      ),
      error = function(e) NULL
    )

    if (is.null(response)) {
      stop("Failed to retrieve the KEGG organism table from /list/organism.")
    }

    organism_df <- utils::read.delim(
      text = response,
      header = FALSE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE
    )

    if (ncol(organism_df) < 4) {
      stop("Unexpected KEGG organism table format returned by /list/organism.")
    }

    #==============================================================
    # Step 4 - Build internal organism table
    #==============================================================
    organism_df <- organism_df[, seq_len(4)]
    colnames(organism_df) <- c("taxon_id", "organism", "species", "phylogeny")
    rownames(organism_df) <- NULL

    ## Separate "Homo sapiens (human)" into species = "Homo sapiens"
    ## and common_name = "human"
    organism_df$common_name <- NA_character_

    has_parentheses <- grepl("\\(.*\\)", organism_df$species)

    organism_df$common_name[has_parentheses] <- trimws(
      sub("^.*\\((.*)\\).*$", "\\1", organism_df$species[has_parentheses])
    )

    organism_df$species[has_parentheses] <- trimws(
      sub("\\s*\\(.*\\)$", "", organism_df$species[has_parentheses])
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
      if (length(query) > 10 && i > 1 && (i - 1) %% 10 == 0) {
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

      if (length(ko_data) == 0) {
        ko_org_list[[ko]] <- character(0)
        next
      }

      ## Extract organism codes from gene IDs (e.g. "eco:b0001" -> "eco")
      gene_ids  <- as.character(ko_data)
      org_codes <- sub(":.*$", "", gene_ids)
      ko_org_list[[ko]] <- unique(org_codes[org_codes != ""])
    }

    ## Report retrieval failures
    if (length(failed_kos) > 0) {
      message(
        "KEGG query failed for ", length(failed_kos), " KO(s): ",
        paste(failed_kos, collapse = ", "),
        "\nThese were skipped."
      )
    }

    #==============================================================
    # Step 6 - Identify KOs with zero organism hits
    #==============================================================
    empty_kos <- names(ko_org_list)[vapply(ko_org_list, length, integer(1)) == 0]

    if (length(empty_kos) > 0) {
      message(
        "Note: ", length(empty_kos), " KO(s) returned no organism links:\n  ",
        paste(empty_kos, collapse = ", ")
      )
    }

    #==============================================================
    # Step 7 - Collect all unique candidate organism codes
    #==============================================================
    all_org_codes <- unique(unlist(ko_org_list, use.names = FALSE))

    if (length(all_org_codes) == 0) {
      warning("No organism candidates could be retrieved for any of the provided KO IDs.")
      return(data.frame(
        organism_name = character(0),
        organism_code = character(0),
        overlap_count = integer(0),
        overlap_hits  = character(0),
        stringsAsFactors = FALSE
      ))
    }

    #==============================================================
    # Step 8 - Build overlap summary (which KOs map to each organism)
    #==============================================================
    out_df <- do.call(rbind, lapply(all_org_codes, function(org_code) {
      matched_kos <- names(ko_org_list)[vapply(
        ko_org_list,
        function(x) org_code %in% x,
        logical(1)
      )]

      data.frame(
        organism_code = org_code,
        overlap_count = length(matched_kos),
        overlap_hits  = paste(matched_kos, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))

    #==============================================================
    # Step 9 - Attach species names and retain prokaryotes only
    #==============================================================
    match_idx <- match(out_df$organism_code, organism_df$organism)

    out_df$organism_name <- ifelse(
      !is.na(match_idx) &
        !is.na(organism_df$species[match_idx]) &
        organism_df$species[match_idx] != "",
      organism_df$species[match_idx],
      out_df$organism_code
    )

    out_df$phylogeny <- organism_df$phylogeny[match_idx]

    ## Filter to prokaryotes only and "hsa"
    out_df <- out_df[
      !is.na(out_df$phylogeny) &
        (grepl("^Prokaryotes;", out_df$phylogeny) | out_df$organism_code == "hsa"),
      ,
      drop = FALSE
    ]

    #==============================================================
    # Step 10 - Final formatting and return
    #==============================================================
    out_df <- out_df[, c("organism_name", "organism_code", "overlap_count", "overlap_hits")]
    out_df <- out_df[order(-out_df$overlap_count, out_df$organism_code), ]
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
    if (length(query) > 10 && i > 1 && (i - 1) %% 10 == 0) {
      Sys.sleep(1)
    }

    ## Retrieve compound entry
    compound_info <- tryCatch(
      KEGGREST::keggGet(cpd),
      error = function(e) NULL
    )

    if (is.null(compound_info) || length(compound_info) == 0 || is.null(compound_info[[1]])) {
      failed_cpds <- c(failed_cpds, cpd)
      next
    }

    ## Extract compound name (first synonym, strip trailing semicolons)
    compound_name <- cpd
    if (!is.null(compound_info[[1]]$NAME) && length(compound_info[[1]]$NAME) > 0) {
      compound_name <- sub(";.*$", "", compound_info[[1]]$NAME[1])
    }

    ## Extract linked reaction IDs
    if (is.null(compound_info[[1]]$REACTION) || length(compound_info[[1]]$REACTION) == 0) {
      no_reaction <- c(no_reaction, cpd)
      next
    }

    reaction_ids <- unique(trimws(unlist(strsplit(compound_info[[1]]$REACTION, " "))))
    reaction_ids <- reaction_ids[reaction_ids != ""]

    if (length(reaction_ids) == 0) {
      no_reaction <- c(no_reaction, cpd)
      next
    }

    ## For each reaction, collect orthology (KO) entries
    for (j in seq_along(reaction_ids)) {
      rid <- reaction_ids[j]

      ## Pause every 10 reaction queries to reduce KEGG API pressure
      if (length(reaction_ids) > 10 && j > 1 && (j - 1) %% 10 == 0) {
        Sys.sleep(1)
      }

      reaction_info <- tryCatch(
        KEGGREST::keggGet(rid),
        error = function(e) NULL
      )

      ## Skip if reaction retrieval fails or has no ORTHOLOGY field
      if (is.null(reaction_info) || length(reaction_info) == 0 || is.null(reaction_info[[1]])) next
      if (is.null(reaction_info[[1]]$ORTHOLOGY) || length(reaction_info[[1]]$ORTHOLOGY) == 0) next

      ko_ids   <- names(reaction_info[[1]]$ORTHOLOGY)
      ko_names <- as.character(reaction_info[[1]]$ORTHOLOGY)
      if (length(ko_ids) == 0) next

      out_rows[[row_index]] <- data.frame(
        compound_id   = cpd,
        compound_name = compound_name,
        KO_ID         = ko_ids,
        KO_name       = ko_names,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  #==============================================================
  # Step 5 - Report retrieval failures and unmatched compounds
  #==============================================================
  if (length(failed_cpds) > 0) {
    message(
      "KEGG query failed for ", length(failed_cpds), " compound(s): ",
      paste(failed_cpds, collapse = ", "),
      "\nThese were skipped."
    )
  }

  if (length(no_reaction) > 0) {
    message(
      "Note: ", length(no_reaction), " compound(s) had no linked reactions:\n  ",
      paste(no_reaction, collapse = ", ")
    )
  }

  ## Build output table before checking unmatched compounds
  if (length(out_rows) == 0) {
    out_df <- data.frame(
      compound_id   = character(0),
      compound_name = character(0),
      KO_ID         = character(0),
      KO_name       = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    out_df <- do.call(rbind, out_rows)
  }

  ## Compounds in input but absent from any output row
  matched_cpds   <- unique(out_df$compound_id)
  unmatched_cpds <- setdiff(query, matched_cpds)

  if (length(unmatched_cpds) > 0) {
    message(
      "Summary: ", length(unmatched_cpds), " of ", length(query),
      " compound(s) yielded no candidate KOs:\n  ",
      paste(unmatched_cpds, collapse = ", ")
    )
  }

  #==============================================================
  # Step 6 - Deduplicate and return
  #==============================================================
  if (nrow(out_df) == 0) {
    warning("No candidate KOs could be retrieved for any of the provided compound IDs.")
    return(out_df)
  }

  out_df <- unique(out_df)
  rownames(out_df) <- NULL

  out_df
}
