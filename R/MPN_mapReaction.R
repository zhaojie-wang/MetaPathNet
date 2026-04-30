#' Run a SPARQL query via MetaNetX endpoint using httr/jsonlite
#'
#' Returns the parsed SPARQL result table as a data.frame.
#'
#' @param endpoint SPARQL endpoint URL.
#' @param query Character string containing the SPARQL query.
#'
#' @return A data.frame with one column per SPARQL variable. If no rows are
#' returned, an empty data.frame is returned.
#'
#' @keywords internal
#' @noRd
query_sparql_httr <- function(endpoint, query) {
  res <- httr::POST(
    url = endpoint,
    body = list(
      query  = query,
      format = "application/sparql-results+json"
    ),
    encode = "form"
  )

  httr::stop_for_status(res)

  ## httr cannot auto-parse application/sparql-results+json
  res_txt <- httr::content(res, as = "text", encoding = "UTF-8")
  res_parsed <- jsonlite::fromJSON(res_txt, simplifyVector = FALSE)

  bindings <- res_parsed$results$bindings

  if (is.null(bindings) || length(bindings) == 0) {
    return(data.frame())
  }

  ## Extract all variable names present in the result head
  vars <- res_parsed$head$vars
  if (is.null(vars) || length(vars) == 0) {
    return(data.frame())
  }

  ## Build a plain data.frame (one column per SPARQL variable)
  out <- do.call(
    rbind,
    lapply(bindings, function(b) {
      row <- lapply(vars, function(v) {
        if (!is.null(b[[v]]) && !is.null(b[[v]]$value)) b[[v]]$value else NA_character_
      })
      row <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
      colnames(row) <- vars
      row
    })
  )

  out
}

#' Convert KEGG reactions to MetaPathNet-style network edges
#'
#' Internal helper that retrieves KEGG reaction entries and converts them into
#' directional edges between compounds and KO terms.
#'
#' @param reaction_ids Character vector of KEGG reaction IDs
#'   (e.g. \code{"R01302"} or \code{"rn:R01302"}).
#'
#' @return
#' A character matrix with three columns:
#' \describe{
#'   \item{source}{Source node ID (compound \code{"cpd:Cxxxxx"} or KO \code{"Kxxxxx"}).}
#'   \item{target}{Target node ID.}
#'   \item{interaction_type}{Edge type, e.g. \code{"k_compound:reversible"}
#'         or \code{"k_compound:irreversible"}.}
#' }
#'
#' @details
#' This function is intended for internal use by \code{MPN_mapReaction()} and
#' assumes valid KEGG reaction identifiers that can be resolved via
#' \pkg{KEGGREST}. It parses the \code{EQUATION} field to infer directionality
#' and extract substrate and product compounds.
#'
#' KEGG Orthology (KO) identifiers are retrieved preferentially from the
#' reaction \code{ORTHOLOGY} field. If no KO is available there, the function
#' falls back to KEGG reaction-to-enzyme and enzyme-to-KO mappings. When KO
#' information cannot be recovered, the reaction is represented as direct
#' compound-to-compound edges.
#'
#' @keywords internal
#' @noRd
reaction_to_matrix <- function(reaction_ids) {
  # 1) Basic checks and ID normalisation
  if (missing(reaction_ids) || length(reaction_ids) == 0) {
    stop("reaction_ids must be a non-empty character vector of KEGG reaction IDs.")
  }

  reaction_ids <- unique(trimws(as.character(reaction_ids)))
  reaction_ids <- reaction_ids[reaction_ids != ""]

  if (length(reaction_ids) == 0) {
    stop("All provided reaction_ids are empty after trimming.")
  }

  # Accept both "R01302" and "rn:R01302", but normalise to "R01302"
  reaction_ids_norm <- gsub("^rn:", "", reaction_ids)

  valid_pattern <- grepl("^R\\d{5}$", reaction_ids_norm)
  if (!all(valid_pattern)) {
    bad <- reaction_ids[!valid_pattern]
    stop(
      "Invalid KEGG reaction ID format detected. IDs must look like 'R01302'. ",
      "First invalid IDs: ", paste(head(bad, 5), collapse = ", "),
      " (total invalid: ", length(bad), ")"
    )
  }

  # 2) Loop over reactions and collect edges
  all_edges <- list()
  edge_index <- 1L

  for (rid in reaction_ids_norm) {
    rid_full <- paste0("rn:", rid)  # KEGGREST prefers "rn:Rxxxxx"

    # Fetch reaction entry from KEGG
    kdat <- tryCatch(
      KEGGREST::keggGet(rid_full),
      error = function(e) NULL
    )
    if (is.null(kdat) || length(kdat) == 0 || is.null(kdat[[1]])) {
      message("Skipping ", rid, ": unable to retrieve KEGG entry.")
      next
    }
    entry <- kdat[[1]]

    if (is.null(entry$EQUATION) || length(entry$EQUATION) == 0) {
      message("Skipping ", rid, ": no EQUATION field in KEGG entry.")
      next
    }

    eq <- entry$EQUATION[1]

    # 3) Parse equation into left/right compounds and direction
    if (grepl("<=>", eq, fixed = TRUE)) {
      direction <- "reversible"
      split_sym <- "<=>"
    } else if (grepl("=>", eq, fixed = TRUE)) {
      direction <- "irreversible"
      split_sym <- "=>"
    } else if (grepl("<=", eq, fixed = TRUE)) {
      direction <- "irreversible"
      split_sym <- "<="
    } else {
      message("Skipping ", rid, ": could not detect direction symbol in EQUATION.")
      next
    }

    parts <- strsplit(eq, split_sym, fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      message("Skipping ", rid, ": unexpected EQUATION format.")
      next
    }

    left_str  <- trimws(parts[1])
    right_str <- trimws(parts[2])

    # Extract KEGG compound IDs
    left_c  <- unique(regmatches(left_str,  gregexpr("C\\d{5}", left_str))[[1]])
    right_c <- unique(regmatches(right_str, gregexpr("C\\d{5}", right_str))[[1]])

    if (length(left_c) == 0 || length(right_c) == 0) {
      message("Skipping ", rid, ": could not detect compounds on both sides of EQUATION.")
      next
    }

    left_cpd  <- paste0("cpd:", left_c)
    right_cpd <- paste0("cpd:", right_c)

    # 4) Map reaction -> KO directly from KEGG reaction entry
    ko_ids <- character(0)

    if (!is.null(entry$ORTHOLOGY) && length(entry$ORTHOLOGY) > 0) {
      ko_candidates <- names(entry$ORTHOLOGY)
      ko_candidates <- unique(ko_candidates[grepl("^K\\d{5}$", ko_candidates)])

      if (length(ko_candidates) > 0) {
        # Keep the first KO only to avoid inflating reaction-derived edges
        ko_ids <- ko_candidates[1]
      }
    }

    if (length(ko_ids) == 0) {
      # 4b) Map reaction -> enzyme -> KO (organism-independent)
      enz_links <- tryCatch(
        KEGGREST::keggLink("enzyme", rid_full),
        error = function(e) NULL
      )
      if (!is.null(enz_links) && length(enz_links) > 0) {
        enz_ids <- unique(as.character(enz_links))

        ko_links <- tryCatch(
          KEGGREST::keggLink("ko", enz_ids),
          error = function(e) NULL
        )
        if (!is.null(ko_links) && length(ko_links) > 0) {
          ko_ids <- unique(sub("^ko:", "", as.character(ko_links)))
        }
      }
    }

    # 5) Build MetaPathNet-style edges for this reaction
    itype <- paste0("k_compound:", direction)

    if (length(ko_ids) > 0) {
      # Edges via KO(s)
      for (ko in ko_ids) {
        # substrates -> KO
        for (cpdL in left_cpd) {
          all_edges[[edge_index]] <- data.frame(
            source           = cpdL,
            target           = ko,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source           = ko,
              target           = cpdL,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }

        # KO -> products
        for (cpdR in right_cpd) {
          all_edges[[edge_index]] <- data.frame(
            source           = ko,
            target           = cpdR,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source           = cpdR,
              target           = ko,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }
      }
    } else {
      # Fallback: compound -> compound edges only
      for (cpdL in left_cpd) {
        for (cpdR in right_cpd) {
          all_edges[[edge_index]] <- data.frame(
            source           = cpdL,
            target           = cpdR,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source           = cpdR,
              target           = cpdL,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }
      }
    }
  }

  # 6) Merge all reactions and clean up
  if (length(all_edges) == 0) {
    warning("No edges could be generated for the provided reaction_ids.")
    out <- matrix(character(0), ncol = 3)
    colnames(out) <- c("source", "target", "interaction_type")
    return(out)
  }

  edge_df <- unique(do.call(rbind, all_edges))

  out <- as.matrix(edge_df[, c("source", "target", "interaction_type")])
  colnames(out) <- c("source", "target", "interaction_type")
  rownames(out) <- NULL

  out
}

#' Convert MetaNetX reaction info into network edges
#'
#' Internal helper that fetches compound-level equations from MetaNetX for
#' unmapped external reactions and returns compound-compound edges.
#'
#' @param mnxr_ids Character vector of MetaNetX reaction IDs
#'   (e.g. \code{"MNXR146590"}).
#'
#' @return
#' A character matrix with three columns:
#' \describe{
#'   \item{source}{Source node ID (compound name from MetaNetX).}
#'   \item{target}{Target node ID.}
#'   \item{interaction_type}{Edge type, here \code{"external"}.}
#' }
#'
#' @details
#' This helper is used when MetaNetX reactions cannot be mapped back to KEGG.
#' It builds edges directly from the MetaNetX stoichiometric equation, using
#' human-readable compound names as nodes and marking edges as
#' \code{"external"} to distinguish them from KEGG-derived interactions.
#'
#' @keywords internal
#' @noRd
reactionInfo_to_matrix <- function(mnxr_ids) {
  # 1) Normalise input
  if (missing(mnxr_ids) || length(mnxr_ids) == 0) {
    stop("mnxr_ids must be a non-empty character vector of MetaNetX reaction IDs (e.g. 'MNXR146590').")
  }
  mnxr_ids <- unique(trimws(as.character(mnxr_ids)))
  mnxr_ids <- mnxr_ids[mnxr_ids != ""]
  if (length(mnxr_ids) == 0) {
    stop("All provided mnxr_ids are empty after trimming.")
  }

  endpoint <- "https://rdf.metanetx.org/sparql"
  all_edges <- list()

  # 2) Loop over MNXR IDs and fetch reaction equation
  for (mnxr_id in mnxr_ids) {
    mnx_reac_uri <- paste0("<https://rdf.metanetx.org/reac/", mnxr_id, ">")

    query_eq <- paste0(
      "PREFIX mnx:  <https://rdf.metanetx.org/schema/> ",
      "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> ",
      "SELECT ?chem ?chem_name ?comp ?comp_name ?coef WHERE { ",
      "  ", mnx_reac_uri, " ?side ?part . ",
      "  ?part mnx:chem ?chem ; ",
      "        mnx:comp ?comp ; ",
      "        mnx:coef ?c . ",
      "  ?chem rdfs:comment ?chem_name . ",
      "  ?comp rdfs:comment ?comp_name . ",
      "  FILTER( ?side IN ( mnx:left , mnx:right )) ",
      "  BIND( IF( ?side = mnx:left, - ?c, ?c ) AS ?coef ) ",
      "}"
    )

    df_eq <- tryCatch(
      query_sparql_httr(endpoint = endpoint, query = query_eq),
      error = function(e) NULL
    )

    if (is.null(df_eq) || nrow(df_eq) == 0) {
      warning("No equation found in MetaNetX for MNXR ID: ", mnxr_id)
      next
    }

    # 3) Separate substrates and products via coefficient sign
    df_eq$coef_num <- as.numeric(df_eq$coef)
    subs  <- df_eq[df_eq$coef_num < 0, , drop = FALSE]
    prods <- df_eq[df_eq$coef_num > 0, , drop = FALSE]

    if (nrow(subs) == 0 || nrow(prods) == 0) {
      warning("Could not clearly separate substrates and products for MNXR ID: ", mnxr_id)
      next
    }

    # 4) Use human-readable compound names as node IDs
    sub_names  <- unique(subs$chem_name)
    prod_names <- unique(prods$chem_name)

    # 5) Build MetaPathNet-style edges: all substrate to product pairs
    edge_df <- expand.grid(
      source = sub_names,
      target = prod_names,
      stringsAsFactors = FALSE
    )
    edge_df$interaction_type <- "external"

    edge_mat <- as.matrix(edge_df[, c("source", "target", "interaction_type")])
    rownames(edge_mat) <- NULL
    all_edges[[mnxr_id]] <- edge_mat
  }

  # 6) Combine all reactions into one matrix
  all_edges <- Filter(Negate(is.null), all_edges)

  if (!length(all_edges)) {
    warning("No valid MetaNetX reactions could be converted.")
    return(matrix(character(0), ncol = 3,
                  dimnames = list(NULL, c("source", "target", "interaction_type"))))
  }

  out <- do.call(rbind, all_edges)
  colnames(out) <- c("source", "target", "interaction_type")
  rownames(out) <- NULL
  out
}

#' Map reaction IDs to MetaPathNet-style networks
#'
#' \code{MPN_mapReaction()} converts reaction identifiers into the standard
#' MetaPathNet 3-column edge list format
#' (\code{source}, \code{target}, \code{interaction_type}).
#'
#' The function supports direct KEGG reaction IDs, as well as external reaction
#' identifiers from Rhea and MetaCyc. For external sources, it uses the
#' MetaNetX SPARQL endpoint as a neutral layer to translate reaction IDs into
#' KEGG reaction IDs where possible. KEGG-mapped reactions are then converted
#' into MetaPathNet-style KO–compound edges. If no KEGG mapping is available,
#' MetaNetX reaction equations are used to generate fallback compound-to-compound
#' edges.
#'
#' @param reaction_ids Character vector of reaction identifiers.
#'
#' @param source Character string indicating the origin of the reaction IDs.
#'   Supported values are:
#'   \itemize{
#'     \item \code{"rhea"}    - Rhea reaction IDs (digits only, no prefix).
#'     \item \code{"metacyc"} - MetaCyc reaction IDs.
#'     \item \code{"kegg"}    - KEGG reaction IDs (e.g. \code{"R01302"} or \code{"rn:R01302"}).
#'   }
#'
#' @details
#' For \code{source = "kegg"}, reaction IDs are converted directly using KEGG
#' reaction entries. Reaction equations are parsed to infer substrates,
#' products, and directionality. KEGG Orthology (KO) identifiers are retrieved
#' preferentially from the reaction \code{ORTHOLOGY} field, and if unavailable,
#' the function falls back to reaction-to-enzyme and enzyme-to-KO mapping.
#'
#' For \code{source = "rhea"} and \code{source = "metacyc"}, the function first
#' queries MetaNetX to identify corresponding MetaNetX reactions, then attempts
#' to recover KEGG reaction IDs. KEGG-mapped reactions follow the same semantics
#' as the metabolic layer of MetaPathNet: compounds are coded as
#' \code{"cpd:Cxxxxx"}, KOs as \code{"Kxxxxx"}, and directionality
#' (\code{reversible} vs \code{irreversible}) is derived from the KEGG
#' \code{EQUATION} field.
#'
#' Reactions without KEGG mapping are represented as compound-to-compound
#' relations, using human-readable compound names as node IDs and
#' \code{"external"} as interaction type. This preserves reaction information
#' while clearly marking the non-KEGG origin of these edges.
#'
#' @return
#' A character matrix with three columns:
#' \describe{
#'   \item{source}{Source node ID (compound, KO, or compound name for
#'         non-KEGG reactions).}
#'   \item{target}{Target node ID.}
#'   \item{interaction_type}{For KEGG-mapped reactions, values such as
#'         \code{"k_compound:reversible"} or
#'         \code{"k_compound:irreversible"}; for reactions without KEGG
#'         mapping, \code{"external"}.}
#' }
#' If no valid reactions can be converted, an empty 3-column matrix is returned
#' and a warning is issued.
#'
#' @examples
#' \donttest{
#' ## NOTE: The following examples require live access to KEGG REST and/or the
#' ## MetaNetX SPARQL endpoint.
#'
#' ## 1) Map a MetaCyc reaction and inspect the returned MetaPathNet-style edges
#' trp_metacyc_edges <- MPN_mapReaction(
#'   reaction_ids = "TRYPTOPHAN-RXN",
#'   source       = "metacyc"
#' )
#'
#' ## 2) Map a Rhea reaction
#' rhea16505_edges <- MPN_mapReaction(
#'   reaction_ids = "16505",
#'   source       = "rhea"
#' )
#'
#' ## 3) Map KEGG reaction IDs directly
#' tma_edges <- MPN_mapReaction(
#'   reaction_ids = c("R05623", "R10285"),
#'   source       = "kegg"
#' )
#' }
#'
#' @seealso
#' \code{\link{MPN_keggNetwork}}, \code{\link{MPN_crossSpeciesNetwork}}
#'
#' @importFrom KEGGREST keggGet keggLink
#' @export
MPN_mapReaction <- function(
    reaction_ids,
    source = c("rhea", "metacyc", "kegg")
) {
  #==============================================================
  # Step 1 - Basic checks, normalisation, and package guards
  #==============================================================
  # Ensure user provided at least one reaction ID
  if (missing(reaction_ids) || length(reaction_ids) == 0) {
    stop("reaction_ids must be a non-empty character vector of reaction IDs.")
  }

  # Clean up IDs: character, unique, trim whitespace, drop empties
  reaction_ids <- unique(trimws(as.character(reaction_ids)))
  reaction_ids <- reaction_ids[reaction_ids != ""]
  if (length(reaction_ids) == 0) {
    stop("All provided reaction_ids are empty after trimming.")
  }

  # Restrict source to supported options
  source <- match.arg(source)

  #==============================================================
  # Step 2 - Direct KEGG reaction mapping
  #==============================================================
  # If KEGG reaction IDs are supplied directly, bypass MetaNetX
  if (source == "kegg") {
    out <- reaction_to_matrix(reaction_ids)

    if (is.null(out) || nrow(out) == 0) {
      warning("No edges could be generated from the provided KEGG reaction IDs; returning empty matrix.")
      out_empty <- matrix(character(0), ncol = 3)
      colnames(out_empty) <- c("source", "target", "interaction_type")
      return(out_empty)
    }

    out <- as.matrix(out)
    colnames(out) <- c("source", "target", "interaction_type")
    rownames(out) <- NULL
    return(out)
  }

  endpoint <- "https://rdf.metanetx.org/sparql"

  # Normalize IDs per source:
  # - Rhea: keep numeric core (e.g. 'rhea:16505' -> '16505')
  # - MetaCyc: use reaction ID directly (e.g. 'TRYPTOPHAN-RXN')
  if (source == "rhea") {
    norm_ids <- sub(".*?(\\d+)$", "\\1", reaction_ids)
  } else {
    norm_ids <- reaction_ids
  }

  #==============================================================
  # Step 3 - Map external reaction IDs -> MetaNetX reaction (MNXR)
  #==============================================================
  # One MNXR label per external reaction (NA if mapping fails)
  mnxr_vec <- rep(NA_character_, length(norm_ids))
  names(mnxr_vec) <- reaction_ids

  # Query MetaNetX for each external reaction and extract MNXR label
  for (i in seq_along(norm_ids)) {
    ext_id  <- reaction_ids[i]
    norm_id <- norm_ids[i]

    # Build the external to MetaNetX SPARQL query
    if (source == "rhea") {
      query_ext_to_mnx <- paste0(
        "PREFIX mnx: <https://rdf.metanetx.org/schema/> ",
        "PREFIX rh:  <http://rdf.rhea-db.org/> ",
        "SELECT DISTINCT ?reaction WHERE { ",
        "  ?reaction a mnx:REAC . ",
        "  ?reaction mnx:reacXref rh:", norm_id, " . ",
        "}"
      )
    } else {  # "metacyc"
      query_ext_to_mnx <- paste0(
        "PREFIX mnx:      <https://rdf.metanetx.org/schema/> ",
        "PREFIX metacycR: <https://identifiers.org/metacyc.reaction:> ",
        "SELECT DISTINCT ?reaction WHERE { ",
        "  ?reaction a mnx:REAC . ",
        "  ?reaction mnx:reacXref metacycR:", norm_id, " . ",
        "}"
      )
    }

    # Run SPARQL query safely
    df_ext <- tryCatch(
      query_sparql_httr(endpoint = endpoint, query = query_ext_to_mnx),
      error = function(e) NULL
    )

    if (is.null(df_ext)) {
      warning("SPARQL query failed for ", source, " ID '", ext_id, "'.")
      next
    }

    if (nrow(df_ext) == 0L) {
      next
    }

    # Use the first mapped MetaNetX reaction URI and extract MNXR ID
    mnx_uri <- df_ext$reaction[1]
    mnxr_id <- sub(".*(MNXR[0-9]+).*", "\\1", mnx_uri)

    if (is.na(mnxr_id) || !nzchar(mnxr_id)) {
      warning("Could not parse MNXR ID for ", source, " ID: ", ext_id)
      next
    }

    mnxr_vec[ext_id] <- mnxr_id
  }

  ext_with_mnxr <- names(mnxr_vec)[!is.na(mnxr_vec)]
  mnxr_all      <- mnxr_vec[ext_with_mnxr]

  # If nothing mapped to MetaNetX at all, abort early with empty matrix
  if (!length(mnxr_all)) {
    warning("No external reactions could be mapped to MetaNetX; returning empty matrix.")
    out_empty <- matrix(character(0), ncol = 3)
    colnames(out_empty) <- c("source", "target", "interaction_type")
    return(out_empty)
  }

  #==============================================================
  # Step 4 - Map MetaNetX reactions (MNXR) -> KEGG reaction IDs
  #==============================================================
  # For each MNXR ID, try to retrieve a KEGG reaction cross-reference
  kegg_vec <- rep(NA_character_, length(mnxr_all))
  names(kegg_vec) <- mnxr_all

  for (mnxr_id in mnxr_all) {
    query_mnxr_to_kegg <- paste0(
      "PREFIX mnx:  <https://rdf.metanetx.org/schema/> ",
      "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> ",
      "PREFIX keggR: <https://identifiers.org/kegg.reaction:> ",
      "SELECT DISTINCT ?kegg_reaction WHERE { ",
      "  ?reac a mnx:REAC . ",
      "  ?reac rdfs:label \"", mnxr_id, "\" . ",
      "  ?reac mnx:reacXref ?kegg_reaction . ",
      "  FILTER(STRSTARTS(STR(?kegg_reaction), STR(keggR:))) ",
      "}"
    )

    df_kegg <- tryCatch(
      query_sparql_httr(endpoint = endpoint, query = query_mnxr_to_kegg),
      error = function(e) NULL
    )

    if (is.null(df_kegg)) {
      warning("SPARQL query MNXR -> KEGG failed for ", mnxr_id, ".")
      next
    }

    if (nrow(df_kegg) == 0L) {
      next
    }

    # Take first KEGG hit and extract reaction ID (Rxxxxx)
    kegg_uri <- df_kegg$kegg_reaction[1]
    kegg_id  <- sub("^.*kegg\\.reaction:([^>]+)>?.*$", "\\1", kegg_uri)

    if (!is.na(kegg_id) && nzchar(kegg_id)) {
      kegg_vec[mnxr_id] <- kegg_id
    }
  }

  #==============================================================
  # Step 5 - Partition MNXR with KEGG vs MNXR without KEGG
  #==============================================================
  # KEGG-mapped reactions
  kegg_ids_unique <- unique(kegg_vec[!is.na(kegg_vec)])

  # MNXR IDs that must be resolved purely from MetaNetX equations
  mnxr_without_kegg <- names(kegg_vec)[is.na(kegg_vec)]

  #==============================================================
  # Step 6 - Build network from KEGG reactions (reaction_to_matrix)
  #==============================================================
  # Use KEGG reaction IDs wherever available to obtain KO-compound edges
  net_kegg <- NULL
  if (length(kegg_ids_unique) > 0L) {
    net_kegg <- reaction_to_matrix(kegg_ids_unique)
  }

  #==============================================================
  # Step 7 - Build network from MetaNetX equations (reactionInfo_to_matrix)
  #==============================================================
  # For MNXR without KEGG mapping, derive compound-compound edges from MNX equation
  net_mnxr <- NULL
  if (length(mnxr_without_kegg) > 0L) {
    net_mnxr <- reactionInfo_to_matrix(mnxr_without_kegg)
  }

  #==============================================================
  # Step 8 - Merge KEGG-derived and MNXR-derived networks and return
  #==============================================================
  # If nothing produced any edges, return an empty MetaPathNet-style matrix
  if (is.null(net_kegg) && is.null(net_mnxr)) {
    warning("No edges could be generated from external reactions; returning empty matrix.")
    out_empty <- matrix(character(0), ncol = 3)
    colnames(out_empty) <- c("source", "target", "interaction_type")
    return(out_empty)
  }

  # Merge both sources when both are non-empty; otherwise keep the one that exists
  if (!is.null(net_kegg) && !is.null(net_mnxr) &&
      nrow(net_kegg) > 0 && nrow(net_mnxr) > 0) {
    out <- unique(rbind(net_kegg, net_mnxr))
  } else if (!is.null(net_kegg) && nrow(net_kegg) > 0) {
    out <- unique(net_kegg)
  } else {
    out <- unique(net_mnxr)
  }

  # Ensure a clean 3-column character matrix as final output
  out <- as.matrix(out)
  colnames(out) <- c("source", "target", "interaction_type")
  rownames(out) <- NULL

  out
}
