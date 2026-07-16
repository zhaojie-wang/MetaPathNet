#' Find KEGG Identifiers for Compounds, Organisms, and KO Terms
#'
#' Map user-supplied compound, organism, or KO queries to KEGG
#' identifiers.
#'
#' @param KEGG_database Character string specifying the KEGG entry type:
#'   \code{"compound"}, \code{"organism"}, or \code{"ko"}.
#' @param searchBy Character string specifying the query type.
#'   For \code{KEGG_database = "compound"}: one of \code{"name"},
#'   \code{"sid"}, \code{"cid"}, \code{"SMILES"}, \code{"InChI"},
#'   or \code{"InChIKey"}.
#'   For \code{KEGG_database = "organism"}: one of \code{"name"} or
#'   \code{"taxon_id"}, where \code{"taxon_id"} refers to the KEGG genome
#'   T number.
#'   For \code{KEGG_database = "ko"}: one of \code{"name"},
#'   \code{"symbol"}, or \code{"ECnumber"}.
#' @param query Vector of queries to be matched. Non-character input is
#'   coerced to \code{character}.
#'
#' @return A \code{data.frame} with the input values and matched KEGG
#'   identifiers:
#'   \itemize{
#'     \item for \code{KEGG_database = "compound"}: columns
#'           \code{Input}, \code{KEGG_ID};
#'     \item for \code{KEGG_database = "organism"}: columns
#'           \code{Input}, \code{Organism_Code};
#'     \item for \code{KEGG_database = "ko"}: columns
#'           \code{Input}, \code{KO_ID}.
#'   }
#'   If no match is found, the corresponding identifier is \code{NA}.
#'   Some inputs may return multiple rows when multiple matches are found.
#'
#' @details
#' For compounds, name queries are first matched against the live KEGG
#' compound table. Unresolved names and structure-based inputs are queried
#' through PubChem, and supported PubChem identifiers are converted to KEGG
#' compound IDs when possible.
#'
#' For organisms, matching is performed against the live KEGG genome table.
#' Scientific and common names are matched exactly without regard to case.
#'
#' For KO terms, the current KEGG KO table is retrieved through the KEGG
#' REST API. Symbol queries may return multiple KO IDs combined with
#' \code{"/"}.
#'
#' @examples
#' ## 1) Map compounds to KEGG compound IDs
#' MPN_keggFinder(
#'   KEGG_database = "compound",
#'   searchBy      = "name",
#'   query         = c("tryptophan", "indole", "serotonin")
#' )
#'
#' ## 2) Map organisms to KEGG organism codes by name
#' ## This example requires live KEGG genome access.
#' \donttest{
#' MPN_keggFinder(
#'   KEGG_database = "organism",
#'   searchBy      = "name",
#'   query         = c("human", "horse", "dragon", "Escherichia coli W")
#' )
#' }
#'
#' ## 3) Map enzyme EC numbers to KEGG KO identifiers
#' MPN_keggFinder(
#'   KEGG_database = "ko",
#'   searchBy      = "ECnumber",
#'   query         = c("4.1.99.1", "1.13.11.27")
#' )
#'
#' @importFrom utils read.delim
#' @importFrom curl curl_fetch_memory new_handle
#' @importFrom jsonlite fromJSON
#' @importFrom webchem get_cid
#' @importFrom KEGGREST keggList
#' @export
MPN_keggFinder <- function(
    KEGG_database = c("compound", "organism", "ko"),
    searchBy = NULL,
    query = NULL
) {
  #==============================================================
  # Step 0 — Select KEGG namespace and validate query
  #==============================================================
  KEGG_database <- match.arg(KEGG_database)

  if (is.null(query)) {
    stop("You must supply 'query'.")
  }

  query_chr <- trimws(as.character(query))
  query_chr <- query_chr[!is.na(query_chr) & nzchar(query_chr)]

  if (length(query_chr) == 0L) {
    stop("All provided query values are missing or empty after trimming.")
  }

  #==============================================================
  # Step 1 — Compound lookup using live KEGG + PubChem fallback
  #==============================================================
  if (KEGG_database == "compound") {
    ## Validate accepted input types for compound lookup
    valid_types <- c("name", "sid", "cid", "SMILES", "InChI", "InChIKey")
    if (is.null(searchBy) || !(searchBy %in% valid_types)) {
      stop(sprintf(
        "searchBy must be one of: %s",
        paste(valid_types, collapse = ", ")
      ))
    }

    out_df <- data.frame(
      Input   = character(0),
      KEGG_ID = character(0),
      stringsAsFactors = FALSE
    )

    unresolved_query <- query_chr

    ## For name input, first try direct KEGG compound name / synonym matching
    if (searchBy == "name") {

      ## Retrieve the live KEGG compound table safely
      compound_res <- tryCatch(
        curl::curl_fetch_memory(
          "https://rest.kegg.jp/list/compound",
          handle = curl::new_handle()
        ),
        error = function(e) NULL
      )

      compound_df <- NULL

      ## Parse the response only if KEGG returned non-empty content
      if (
        !is.null(compound_res) &&
        isTRUE(compound_res$status_code == 200L) &&
        length(compound_res$content) > 0L
      ) {
        response <- rawToChar(compound_res$content)

        if (nzchar(trimws(response))) {
          compound_df <- tryCatch(
            utils::read.delim(
              text = response,
              header = FALSE,
              sep = "\t",
              quote = "",
              stringsAsFactors = FALSE
            ),
            error = function(e) NULL
          )
        }
      }

      ## Keep only valid KEGG compound tables
      if (!is.null(compound_df) && ncol(compound_df) >= 2L) {
        compound_df <- compound_df[, seq_len(2), drop = FALSE]
        colnames(compound_df) <- c("KEGG_ID", "common_names")
        rownames(compound_df) <- NULL
      } else {
        compound_df <- NULL
      }

      ## Match user names against KEGG compound names and synonyms
      if (!is.null(compound_df)) {
        resolved_idx <- logical(length(query_chr))

        for (i in seq_along(query_chr)) {
          q <- query_chr[i]

          hit_idx <- vapply(compound_df$common_names, function(x) {
            syns <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
            any(tolower(syns) == tolower(q))
          }, logical(1))

          if (any(hit_idx)) {
            tmp <- data.frame(
              Input   = rep(q, sum(hit_idx)),
              KEGG_ID = sub("^cpd:", "", compound_df$KEGG_ID[hit_idx]),
              stringsAsFactors = FALSE
            )

            out_df <- rbind(out_df, tmp)
            resolved_idx[i] <- TRUE
          }
        }

        ## Only unresolved names continue to the PubChem fallback route
        unresolved_query <- query_chr[!resolved_idx]
      }
    }

    ## For direct SID input, convert to KEGG immediately
    if (searchBy == "sid") {
      unresolved_query <- character(0)

      for (q in query_chr) {
        sid_num <- suppressWarnings(as.numeric(q))

        ## Keep only the SID range used for KEGG conversion in this workflow
        if (is.na(sid_num) || sid_num >= 15000) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        conv_res <- tryCatch(
          KEGGREST::keggConv("compound", paste0("pubchem:", q)),
          error = function(e) NULL
        )

        if (is.null(conv_res) || length(conv_res) == 0) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        kegg_ids <- unique(sub("^(cpd:|compound:)", "", unname(conv_res)))
        kegg_ids <- kegg_ids[!is.na(kegg_ids) & kegg_ids != ""]

        if (length(kegg_ids) > 0) {
          out_df <- rbind(
            out_df,
            data.frame(
              Input   = rep(q, length(kegg_ids)),
              KEGG_ID = kegg_ids,
              stringsAsFactors = FALSE
            )
          )
        } else {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
        }
      }
    }

    ## For direct CID input, convert CID -> SID -> KEGG
    if (searchBy == "cid") {
      unresolved_query <- character(0)

      for (q in query_chr) {
        url <- paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
          q,
          "/sids/JSON?sids_type=all"
        )

        res <- tryCatch(
          jsonlite::fromJSON(url, simplifyVector = FALSE),
          error = function(e) NULL
        )

        if (is.null(res) ||
            is.null(res$InformationList) ||
            is.null(res$InformationList$Information) ||
            length(res$InformationList$Information) == 0 ||
            is.null(res$InformationList$Information[[1]]$SID)) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        sid_vec <- as.character(res$InformationList$Information[[1]]$SID)
        sid_num <- suppressWarnings(as.numeric(sid_vec))
        sid_vec <- sid_vec[!is.na(sid_num) & sid_num < 15000]

        if (length(sid_vec) == 0) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        kegg_hits <- character(0)

        for (sid in sid_vec) {
          conv_res <- tryCatch(
            KEGGREST::keggConv("compound", paste0("pubchem:", sid)),
            error = function(e) NULL
          )

          if (is.null(conv_res) || length(conv_res) == 0) {
            next
          }

          kegg_ids <- unique(sub("^(cpd:|compound:)", "", unname(conv_res)))
          kegg_ids <- kegg_ids[!is.na(kegg_ids) & kegg_ids != ""]

          if (length(kegg_ids) > 0) {
            kegg_hits <- c(kegg_hits, kegg_ids)
          }
        }

        kegg_hits <- unique(kegg_hits)

        if (length(kegg_hits) > 0) {
          out_df <- rbind(
            out_df,
            data.frame(
              Input   = rep(q, length(kegg_hits)),
              KEGG_ID = kegg_hits,
              stringsAsFactors = FALSE
            )
          )
        } else {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
        }
      }
    }

    ## PubChem fallback for unresolved names, or all structure-based input types
    if (searchBy %in% c("name", "SMILES", "InChI", "InChIKey")) {
      if (searchBy == "name") {
        from_pubchem <- "name"
      } else if (searchBy == "SMILES") {
        from_pubchem <- "smiles"
      } else if (searchBy == "InChI") {
        from_pubchem <- "inchi"
      } else if (searchBy == "InChIKey") {
        from_pubchem <- "inchikey"
      }

      if (searchBy != "name") {
        unresolved_query <- query_chr
      }

      for (q in unresolved_query) {
        ## Retrieve all PubChem CID candidates for the current input
        cid_res <- tryCatch(
          webchem::get_cid(
            query   = q,
            from    = from_pubchem,
            match   = "all",
            verbose = FALSE
          ),
          error = function(e) NULL
        )

        ## If PubChem lookup fails completely, keep the input with NA
        if (is.null(cid_res)) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        cid_res <- as.data.frame(cid_res, stringsAsFactors = FALSE)

        ## If no CID candidate is found, keep the input with NA
        if (nrow(cid_res) == 0) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        ## Extract CID values as robustly as possible from the returned table
        cid_col <- intersect(c("cid", "CID"), colnames(cid_res))

        if (length(cid_col) == 0) {
          fallback_cols <- setdiff(colnames(cid_res), c("query", "Query", "input", "Input"))
          if (length(fallback_cols) == 0) {
            out_df <- rbind(
              out_df,
              data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
            )
            next
          }
          cid_values <- cid_res[[fallback_cols[1]]]
        } else {
          cid_values <- cid_res[[cid_col[1]]]
        }

        cid_values <- unique(trimws(as.character(cid_values)))
        cid_values <- cid_values[!is.na(cid_values) & cid_values != ""]

        ## If CID extraction fails, keep the input with NA
        if (length(cid_values) == 0) {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
          next
        }

        ## Collect KEGG compound IDs recovered from valid SID candidates
        kegg_hits <- character(0)

        for (cid in cid_values) {
          ## Convert each CID to PubChem SIDs
          url <- paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
            cid,
            "/sids/JSON?sids_type=all"
          )

          res <- tryCatch(
            jsonlite::fromJSON(url, simplifyVector = FALSE),
            error = function(e) NULL
          )

          if (is.null(res) ||
              is.null(res$InformationList) ||
              is.null(res$InformationList$Information) ||
              length(res$InformationList$Information) == 0 ||
              is.null(res$InformationList$Information[[1]]$SID)) {
            next
          }

          sid_vec <- as.character(res$InformationList$Information[[1]]$SID)

          ## Keep only SIDs used for KEGG conversion in this workflow
          sid_num <- suppressWarnings(as.numeric(sid_vec))
          sid_vec <- sid_vec[!is.na(sid_num) & sid_num < 15000]

          if (length(sid_vec) == 0) {
            next
          }

          ## Convert each SID to KEGG compound ID
          for (sid in sid_vec) {
            conv_res <- tryCatch(
              KEGGREST::keggConv("compound", paste0("pubchem:", sid)),
              error = function(e) NULL
            )

            if (is.null(conv_res) || length(conv_res) == 0) {
              next
            }

            kegg_ids <- unique(sub("^(cpd:|compound:)", "", unname(conv_res)))
            kegg_ids <- kegg_ids[!is.na(kegg_ids) & kegg_ids != ""]

            if (length(kegg_ids) > 0) {
              kegg_hits <- c(kegg_hits, kegg_ids)
            }
          }
        }

        kegg_hits <- unique(kegg_hits)

        ## Keep mapped KEGG rows when available; otherwise keep one NA row
        if (length(kegg_hits) > 0) {
          out_df <- rbind(
            out_df,
            data.frame(
              Input   = rep(q, length(kegg_hits)),
              KEGG_ID = kegg_hits,
              stringsAsFactors = FALSE
            )
          )
        } else {
          out_df <- rbind(
            out_df,
            data.frame(Input = q, KEGG_ID = NA_character_, stringsAsFactors = FALSE)
          )
        }
      }
    }

    rownames(out_df) <- NULL
    return(out_df)
  }

  #==============================================================
  # Step 2 — Organism lookup using live KEGG genome table
  #==============================================================
  if (KEGG_database == "organism") {
    ## Validate accepted input types for organism lookup
    valid_types <- c("name", "taxon_id")
    if (is.null(searchBy) || !(tolower(searchBy) %in% valid_types)) {
      stop(sprintf(
        "searchBy must be one of: %s",
        paste(valid_types, collapse = ", ")
      ))
    }

    searchBy <- tolower(searchBy)

    ## Retrieve the live KEGG genome table
    genome <- tryCatch(
      KEGGREST::keggList("genome"),
      error = function(e) NULL
    )

    if (is.null(genome) || length(genome) == 0L || is.null(names(genome))) {
      stop("Failed to retrieve the KEGG genome table.")
    }

    genome_value <- trimws(as.character(genome))

    ## Keep only rows matching the current KEGG genome format
    valid_genome <- !is.na(names(genome)) &
      nzchar(names(genome)) &
      !is.na(genome_value) &
      nzchar(genome_value) &
      grepl(";", genome_value, fixed = TRUE)

    if (!any(valid_genome)) {
      stop("KEGG genome table returned no parseable rows.")
    }

    genome <- genome[valid_genome]
    genome_value <- genome_value[valid_genome]

    ## Build the organism lookup table from the genome response
    organism_df <- data.frame(
      taxon_id = names(genome),
      organism = trimws(sub(";.*$", "", genome_value)),
      species  = trimws(sub("^[^;]+;\\s*", "", genome_value)),
      stringsAsFactors = FALSE
    )

    ## Remove any incomplete rows after parsing
    organism_df <- organism_df[
      nzchar(organism_df$taxon_id) &
        nzchar(organism_df$organism) &
        nzchar(organism_df$species),
      ,
      drop = FALSE
    ]

    if (nrow(organism_df) == 0L) {
      stop("KEGG genome table returned no valid organism entries.")
    }

    rownames(organism_df) <- NULL

    ## Split "Homo sapiens (human)" into scientific and common names
    organism_df$common_name <- NA_character_

    has_parentheses <- grepl("\\(.*\\)", organism_df$species)

    organism_df$common_name[has_parentheses] <- trimws(
      sub("^.*\\((.*)\\).*$", "\\1", organism_df$species[has_parentheses])
    )

    organism_df$species[has_parentheses] <- trimws(
      sub("\\s*\\(.*\\)$", "", organism_df$species[has_parentheses])
    )

    out_df <- data.frame(
      Input         = character(0),
      Organism_Code = character(0),
      stringsAsFactors = FALSE
    )

    for (q in query_chr) {
      if (searchBy == "name") {
        ## Match both scientific name and common name
        idx_species <- which(
          !is.na(organism_df$species) &
            tolower(organism_df$species) == tolower(q)
        )

        idx_common <- which(
          !is.na(organism_df$common_name) &
            tolower(organism_df$common_name) == tolower(q)
        )

        idx <- unique(c(idx_species, idx_common))

      } else if (searchBy == "taxon_id") {
        idx <- which(
          !is.na(organism_df$taxon_id) &
            as.character(organism_df$taxon_id) == as.character(q)
        )
      }

      if (length(idx) == 0) {
        out_df <- rbind(
          out_df,
          data.frame(
            Input = q,
            Organism_Code = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      } else {
        out_df <- rbind(
          out_df,
          data.frame(
            Input = rep(q, length(idx)),
            Organism_Code = organism_df$organism[idx],
            stringsAsFactors = FALSE
          )
        )
      }
    }

    rownames(out_df) <- NULL
    return(out_df)
  }

  #==============================================================
  # Step 3 — KO lookup using live KEGG REST table
  #==============================================================
  if (KEGG_database == "ko") {
    ## Validate accepted input types for KO lookup
    valid_types <- c("name", "ECnumber", "symbol")
    if (is.null(searchBy) || !(searchBy %in% valid_types)) {
      stop(sprintf(
        "searchBy must be one of: %s",
        paste(valid_types, collapse = ", ")
      ))
    }

    ## Retrieve KO descriptions from KEGG REST API
    ko_res <- tryCatch(
      curl::curl_fetch_memory(
        "https://rest.kegg.jp/list/ko",
        handle = curl::new_handle()
      ),
      error = function(e) NULL
    )

    ## Stop before parsing if KEGG returned an empty or failed response
    if (
      is.null(ko_res) ||
      !isTRUE(ko_res$status_code == 200L) ||
      length(ko_res$content) == 0L
    ) {
      stop("Failed to retrieve the KEGG KO table.")
    }

    ko_response <- rawToChar(ko_res$content)

    if (!nzchar(trimws(ko_response))) {
      stop("KEGG KO table returned an empty response.")
    }

    ## Parse the KO table only after confirming non-empty content
    ko_df <- tryCatch(
      utils::read.delim(
        text = ko_response,
        header = FALSE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(ko_df) || ncol(ko_df) < 2L || nrow(ko_df) == 0L) {
      stop("Unexpected KEGG KO table format.")
    }

    ## Keep the standard KO ID and description columns
    ko_df <- ko_df[, seq_len(2), drop = FALSE]
    colnames(ko_df) <- c("KO", "Description")

    ## Pre-compute symbol, name, and EC-number fields from KEGG descriptions
    ko_df$Symbols <- sub("^(.*?)[,; ].*", "\\1", ko_df$Description)
    ko_df$AllSymbols <- lapply(ko_df$Description, function(desc) {
      sub_part <- sub("^([^;]+).*", "\\1", desc)
      unlist(strsplit(sub_part, "[,; ]+"))
    })
    ko_df$Name <- sub("^[^;]+; *([^[]+).*", "\\1", ko_df$Description)
    ko_df$ECs <- regmatches(
      ko_df$Description,
      gregexpr("EC:[0-9.\\-]+", ko_df$Description)
    )

    ## Perform KO lookup for each input query
    out <- vapply(query_chr, function(val) {
      if (is.na(val) || val == "") {
        return(NA_character_)
      }

      if (searchBy == "symbol") {
        idx <- which(vapply(
          ko_df$AllSymbols,
          function(x) any(tolower(x) == tolower(val)),
          logical(1)
        ))

        if (length(idx) == 0) {
          return(NA_character_)
        }

        ko_matches <- unique(ko_df$KO[idx])

        if (length(ko_matches) > 1) {
          message(
            "Multiple KO matches found for '", val, "': ",
            paste(ko_matches, collapse = "/")
          )
        }

        return(paste(ko_matches, collapse = "/"))

      } else if (searchBy == "name") {
        idx <- which(grepl(tolower(val), tolower(ko_df$Name), fixed = TRUE))

      } else if (searchBy == "ECnumber") {
        val_mod <- ifelse(startsWith(val, "EC:"), val, paste0("EC:", val))
        idx <- which(vapply(
          ko_df$ECs,
          function(ec) any(tolower(ec) == tolower(val_mod)),
          logical(1)
        ))
      }

      if (length(idx) > 0) ko_df$KO[idx[1]] else NA_character_
    }, character(1))

    out_df <- data.frame(
      Input = query_chr,
      KO_ID = out,
      stringsAsFactors = FALSE
    )
    rownames(out_df) <- NULL
    return(out_df)
  }
}
