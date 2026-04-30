#' Build KEGG Metabolic and Signaling Network
#'
#' Constructs a directed network (3-column edge list) from user-specified KEGG
#' metabolic and/or signaling pathways. KEGG KGML files are downloaded and
#' parsed, networks are merged, organism-specific gene IDs in signaling
#' pathways are converted to KEGG Orthology (KO) IDs, self-loops are removed,
#' and interaction types are annotated.
#'
#' @param metabo_paths Character vector of KEGG pathway IDs for metabolic
#'   pathways (e.g. \code{"hsa00380"}, \code{"eco00380"}). Use \code{NULL} to
#'   exclude the metabolic layer.
#' @param signaling_paths Character vector of KEGG pathway IDs for signaling
#'   pathways (e.g. \code{"hsa04060"}, \code{"hsa04630"}). Use \code{NULL} to
#'   exclude the signaling layer.
#'
#' @return
#' A \code{data.frame} with 3 columns:
#' \describe{
#'   \item{source}{Source node (compound or KO/gene).}
#'   \item{target}{Target node (KO/gene or compound).}
#'   \item{interaction_type}{Interaction type
#'         (e.g. metabolic reaction, activation, inhibition, phosphorylation).}
#' }
#'
#' @details
#' \itemize{
#'   \item All pathways supplied in \code{metabo_paths} and
#'         \code{signaling_paths} must belong to the same organism
#'         (same prefix, e.g. \code{"hsa"}, \code{"eco"}, \code{"bsu"}).
#'   \item If both metabolic and signaling pathways are provided, the output
#'         network is the union of both layers.
#'   \item Invalid or unavailable pathway IDs are skipped with a message, and
#'         the network is built from valid pathways only.
#' }
#'
#' @seealso
#' \code{\link{MPN_getPathIDs}} for pathway ID retrieval.
#'
#' @examples
#' \dontrun{
#' ## Requires extensive KEGG API querying; runtime depends on server response
#' ## Example 1 — Human tryptophan metabolism + immune signaling
#'
#' ## Metabolic pathway: Tryptophan metabolism (human)
#' metabo_paths_hsa <- c("hsa00380")
#'
#' ## Selected immune / inflammatory signaling pathways (human)
#' signaling_paths_hsa <- c(
#'   "hsa04060",  # Cytokine–cytokine receptor interaction
#'   "hsa04630",  # JAK–STAT signaling pathway
#'   "hsa04064",  # NF-kappa B signaling pathway
#'   "hsa04660",  # T cell receptor signaling pathway
#'   "hsa04659"   # Th17 cell differentiation
#' )
#'
#' net_trp_hsa <- MPN_keggNetwork(
#'   metabo_paths    = metabo_paths_hsa,
#'   signaling_paths = signaling_paths_hsa
#' )
#'
#' head(net_trp_hsa)
#'
#' ## Example 2 — E. coli tryptophan metabolism + biosynthesis
#'
#' metabo_paths_eco <- c(
#'   "eco00380",  # Tryptophan metabolism (E. coli)
#'   "eco00400"   # Phenylalanine, tyrosine and tryptophan biosynthesis (E. coli)
#' )
#'
#' net_trp_eco <- MPN_keggNetwork(
#'   metabo_paths    = metabo_paths_eco,
#'   signaling_paths = NULL
#' )
#'
#' head(net_trp_eco)
#' }
#'
#' @importFrom KEGGREST keggList
#' @importFrom RCurl getURL
#' @importFrom stats aggregate
#' @export
MPN_keggNetwork <- function(
    metabo_paths = NULL,
    signaling_paths = NULL
) {
  # --- Safety: initialize expected global object ---
  interaction_type <- NULL
  lookup <- NULL

  #==============================================================
  # Step 1. Validate input paths
  #==============================================================
  if (is.null(metabo_paths) && is.null(signaling_paths)) {
    stop("At least one of metabo_paths or signaling_paths must be provided")
  }

  metabo_paths    <- unique(metabo_paths)
  signaling_paths <- unique(signaling_paths)
  input_paths     <- c(metabo_paths, signaling_paths)

  #==============================================================
  # Step 2. Detect organism code from KEGG path IDs
  #==============================================================
  first_path <- input_paths[1]
  use4 <- suppressWarnings(is.na(as.numeric(substr(first_path, 4, 4))))
  org_len <- if (use4) 4 else 3
  organism_code <- unique(substr(input_paths, 1, org_len))

  if (length(organism_code) > 1) {
    stop("All paths must belong to the same organism: check path IDs")
  }

  #==============================================================
  # Step 3. Get list of available KEGG pathways for this organism
  #==============================================================
  lines <- tryCatch(
    KEGGREST::keggList("pathway", organism_code),
    error = function(e) NULL
  )
  if (is.null(lines)) {
    stop("Failed to fetch KEGG pathways for organism ", organism_code)
  }
  all_paths <- names(lines)

  #==============================================================
  # Step 4. Validate that input path IDs exist in KEGG
  #==============================================================
  if (!is.null(metabo_paths) && !any(metabo_paths %in% all_paths)) {
    stop("None of the metabo_paths are valid for organism ", organism_code)
  }
  if (!is.null(signaling_paths) && !any(signaling_paths %in% all_paths)) {
    stop("None of the signaling_paths are valid for organism ", organism_code)
  }

  paths_included <- setNames(rep(TRUE, length(input_paths)), input_paths)

  #==============================================================
  # Step 5. Build metabolic network
  #==============================================================
  metabolic_table <- NULL
  if (length(metabo_paths)) {
    message("[metabolic] Building network for ", organism_code)

    metab_list <- lapply(metabo_paths, get_metabonet, all_paths, organism_code)
    names(metab_list) <- metabo_paths

    # Remove empty/invalid results
    empty <- vapply(metab_list, function(x) length(x) <= 1, logical(1))
    if (any(empty)) {
      metab_list <- metab_list[!empty]
      paths_included[metabo_paths[empty]] <- FALSE
    }

    if (length(metab_list)) {
      metabolic_table <- metabolic_matrix(names(metab_list), metab_list, organism_code)

      if (!is.null(metabolic_table)) {
        idx <- grepl("_", metabolic_table)
        metabolic_table[idx] <- substring(metabolic_table[idx], 11)
      } else {
        warning("No valid metabolic network could be built from the provided metabolic paths")
      }
    } else {
      warning("No valid metabolic paths to build network")
    }
  }

  #==============================================================
  # Step 6. Build signaling network
  #==============================================================
  signaling_table <- NULL
  if (length(signaling_paths)) {
    message("[signaling] Building network for ", organism_code)

    network_list <- lapply(signaling_paths, get_signalnet, all_paths)
    names(network_list) <- signaling_paths

    # Remove all invalid: those which are vectors (e.g., "bad_path") or NULL
    is_invalid <- vapply(network_list, function(x) is.null(x) || is.vector(x), logical(1))
    if (any(is_invalid)) {
      network_list <- network_list[!is_invalid]
      paths_included[signaling_paths[is_invalid]] <- FALSE
    }

    if (length(network_list) == 0) {
      warning("No valid signaling networks to build network")
    } else {
      # Merge graphs
      global_network_all <- KEGGgraph::mergeGraphs(network_list)
      signaling_table <- signaling_matrix(global_network_all)

      if (is.null(signaling_table)) {
        warning("No valid signaling edges to build network")
      } else {
        # KO conversion
        url <- paste0("https://rest.kegg.jp/link/ko/", organism_code)
        response_ko <- tryCatch(
          RCurl::getURL(url),
          error = function(e) NULL
        )

        if (is.null(response_ko)) {
          stop("Failed to retrieve KO mappings for organism ", organism_code)
        }

        ko_table <- convertTable(response_ko)

        if (!is.matrix(ko_table) || ncol(ko_table) < 2) {
          stop("Failed to parse KO mappings for organism ", organism_code)
        }

        ko_table[, 2] <- substr(ko_table[, 2], 4, nchar(ko_table[, 2]))
        lookup <- setNames(ko_table[, 2], ko_table[, 1])

        for (i in seq_len(2)) {
          signaling_table[, i] <- vapply(
            signaling_table[, i],
            function(x) {
              if (x %in% names(lookup)) lookup[[x]] else x
            },
            character(1)
          )
        }
      }
    }
  }

  #==============================================================
  # Step 6b. Retrieve signaling interaction types (inline implementation)
  #==============================================================

  interaction_type <- NULL
  if (length(signaling_paths)) {
    message("[signaling] Retrieving interaction types from KEGG KGML")

    valid_paths <- intersect(all_paths, signaling_paths)

    if (length(valid_paths) == 0) {
      warning("No valid signaling paths to fetch interaction types.")
    } else {
      # Fetch and parse KGML relations for all signaling paths
      interaction_list <- lapply(valid_paths, function(path) {
        file <- paste0("https://rest.kegg.jp/get/", path, "/kgml")
        pathway_xml <- tryCatch(
          RCurl::getURL(file),
          error = function(e) NULL
        )
        if (is.null(pathway_xml)) return(NULL)

        df <- tryCatch(
          KEGGgraph::parseKGML2DataFrame(pathway_xml, reactions = FALSE),
          error = function(e) NULL
        )
        if (is.null(df) || nrow(df) == 0) return(NULL)

        cn <- colnames(df)
        src <- if ("entry1" %in% cn) "entry1" else if ("from" %in% cn) "from" else NA
        tgt <- if ("entry2" %in% cn) "entry2" else if ("to" %in% cn) "to" else NA
        sub <- if ("subtype" %in% cn) "subtype" else if ("type" %in% cn) "type" else NA
        if (is.na(src) || is.na(tgt)) return(NULL)

        out <- data.frame(
          source = as.character(df[[src]]),
          target = as.character(df[[tgt]]),
          interaction_type = if (!is.na(sub)) as.character(df[[sub]]) else "unknown",
          stringsAsFactors = FALSE
        )
        out$interaction_type <- trimws(out$interaction_type)
        out
      })

      interaction_list <- Filter(Negate(is.null), interaction_list)

      if (length(interaction_list) == 0) {
        warning("No interaction relations found in KGML for provided signaling paths.")
      } else {
        # Combine all interactions into one dataframe
        interaction_df <- unique(do.call(rbind, interaction_list))
        if (nrow(interaction_df) > 0) {
          # Collapse multiple subtypes for the same source–target pair
          collapsed <- stats::aggregate(
            interaction_type ~ source + target,
            data = interaction_df,
            FUN = function(x) {
              paste(sort(unique(gsub(" ", "-", x))), collapse = "/")
            }
          )
          collapsed$interaction_type <- gsub("compound", "indirect-compound", collapsed$interaction_type, fixed = TRUE)

          # --- Build a robust lookup that matches with or without organism prefix ---
          if (!is.null(lookup)) {
            # lookup: e.g., names "bsu:BSU03440" -> values "K12345"
            lk_pref  <- lookup
            lk_nopfx <- setNames(unname(lk_pref), sub("^[a-z]{3,4}:", "", names(lk_pref)))
            lk <- c(lk_pref, lk_nopfx)

            # Convert species-specific gene IDs → KO IDs (both with or without prefix)
            collapsed$source <- ifelse(collapsed$source %in% names(lk),
                                       unname(lk[collapsed$source]),
                                       collapsed$source)
            collapsed$target <- ifelse(collapsed$target %in% names(lk),
                                       unname(lk[collapsed$target]),
                                       collapsed$target)
          }

          # --- Final clean-up and conversion ---
          collapsed$source <- sub("^ko:", "", collapsed$source)
          collapsed$target <- sub("^ko:", "", collapsed$target)
          collapsed$source <- sub("^[a-z]{3,4}:", "", collapsed$source)
          collapsed$target <- sub("^[a-z]{3,4}:", "", collapsed$target)

          interaction_type <- as.matrix(collapsed[, c("source", "target", "interaction_type")])
          colnames(interaction_type) <- c("source", "target", "interaction_type")
          rownames(interaction_type) <- NULL
        } else {
          warning("Interaction dataframe is empty after combining KGML relations.")
          interaction_type <- NULL
        }
      }
    }
  }

  #==============================================================
  # Step 7. Combine metabolic & signaling networks
  #==============================================================
  if (is.null(metabolic_table) && is.null(signaling_table)) {
    stop("Unable to build any network: both tables are empty")
  }

  if (!is.null(metabolic_table) && !is.null(signaling_table)) {
    net <- unique(rbind(metabolic_table, signaling_table))
  } else if (!is.null(metabolic_table)) {
    net <- unique(metabolic_table)
  } else {
    net <- unique(signaling_table)
  }

  #==============================================================
  # Step 8. Remove self-loops
  #==============================================================
  if (ncol(net) >= 2) {
    bad <- which(net[,1] == net[,2])
    if (length(bad)) {
      net <- net[-bad, , drop = FALSE]
    }
  }

  #==============================================================
  # Step 9. Map interaction types
  #==============================================================
  message("[interactions] Mapping interaction types")
  interactions <- map_interactions(net, interaction_type)
  interactions <- unique(interactions)

  #==============================================================
  # Step 10. Report network stats
  #==============================================================
  network_features(interactions[, seq_len(2)])

  dropped <- names(paths_included)[!paths_included]
  if (length(dropped)) warning("Excluded paths: ", paste(dropped, collapse = ","))

  return(interactions)
}
