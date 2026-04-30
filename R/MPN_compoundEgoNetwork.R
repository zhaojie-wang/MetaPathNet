#' Map compounds to organism-specific pathways
#'
#' Internal helper that identifies organism-specific pathways associated
#' with a set of input compounds. Used by \code{\link{MPN_compoundEgoNetwork}}
#' as a preprocessing step.
#'
#' @param compounds Character vector of compound identifiers.
#' @param organism_codes Character vector of KEGG organism codes.
#'
#' @return
#' A named list where each element corresponds to an organism and contains
#' the identifiers of pathways associated with the input compounds.
#'
#' @keywords internal
#' @noRd
MPN_mapCompoundsToPathways <- function(compounds, organism_codes) {
  #==============================================================
  # Step 0 - Basic input checks
  #==============================================================
  # Validate compound and organism inputs
  if (missing(compounds) || length(compounds) == 0) {
    stop("You must supply at least one compound ID in 'compounds'.")
  }
  if (missing(organism_codes) || length(organism_codes) == 0) {
    stop("You must supply at least one KEGG organism code in 'organism_codes'.")
  }

  compounds <- unique(trimws(as.character(compounds)))
  compounds <- compounds[!is.na(compounds) & compounds != ""]

  organism_codes <- unique(trimws(as.character(organism_codes)))
  organism_codes <- organism_codes[!is.na(organism_codes) & organism_codes != ""]

  #==============================================================
  # Step 1 - Normalise compound IDs and build ID variants
  #==============================================================
  # Extract uppercase KEGG compound core IDs (without prefix)
  core_ids <- toupper(gsub("^(cpd:)?C", "", compounds, ignore.case = TRUE))
  # Create unprefixed Cxxxxx form
  comp_C <- paste0("C", core_ids)
  # Create prefixed cpd:Cxxxxx form
  comp_cpd <- paste0("cpd:", comp_C)
  # Collect all variants for later matching
  all_comp_variants <- unique(c(comp_C, comp_cpd))

  #==============================================================
  # Step 2 - Compound-centric prefilter
  #==============================================================
  # Initialise mapping from compounds to generic pathways
  compound_to_maps <- list()

  for (i in seq_along(comp_C)) {
    cmp <- comp_C[i]
    message("Prefilter: fetching PATHWAY information for compound: ", cmp)

    ## Pause every 10 compounds to reduce KEGG API pressure
    if (i %% 10 == 0) {
      Sys.sleep(1)
    }

    # Retrieve KEGG entry for this compound
    info <- tryCatch(
      KEGGREST::keggGet(paste0("cpd:", cmp)),
      error = function(e) NULL
    )

    # Extract generic pathway IDs (mapXXXXX), if available
    if (is.null(info) || length(info) == 0L || is.null(info[[1]]$PATHWAY)) {
      compound_to_maps[[cmp]] <- character(0)
    } else {
      these_paths <- names(info[[1]]$PATHWAY)
      these_paths <- these_paths[grepl("^map\\d{5}$", these_paths)]
      compound_to_maps[[cmp]] <- unique(these_paths)
    }
  }

  # Union of all candidate map IDs
  all_candidate_maps <- unique(unlist(compound_to_maps))

  #==============================================================
  # Step 3 - Map generic maps to organism-specific pathway IDs
  #==============================================================
  # Initialise list of candidate pathways per organism
  candidate_paths_by_org <- list()

  for (org in organism_codes) {
    message("Organism: ", org, " prefiltering candidate pathways")

    # Retrieve organism-specific pathway list
    pw_vec <- tryCatch(
      KEGGREST::keggList("pathway", org),
      error = function(e) NULL
    )

    if (is.null(pw_vec) || length(pw_vec) == 0L) {
      message("  No pathway list returned for organism: ", org)
      candidate_paths_by_org[[org]] <- character(0)
      next
    }

    # Convert "path:org00010" -> "org00010"
    all_ids   <- names(pw_vec)
    org_paths <- sub("^path:", "", all_ids)

    # Convert organism-specific IDs to generic mapXXXXX
    org_as_map <- sub(paste0("^", org), "map", org_paths)

    # Find intersection with candidate generic pathways
    matched_maps      <- intersect(org_as_map, all_candidate_maps)
    # Convert back to organism-specific IDs
    matched_org_paths <- sub("^map", org, matched_maps)

    candidate_paths_by_org[[org]] <- unique(matched_org_paths)
  }

  #==============================================================
  # Step 4 - Confirmation via KGML and KEGGgraph
  #==============================================================
  # Initialise confirmed pathway list per organism
  confirmed_paths_by_org <- list()

  for (org in names(candidate_paths_by_org)) {
    org_paths <- candidate_paths_by_org[[org]]

    if (length(org_paths) == 0) {
      message("Organism: ", org, " has no candidate pathways after prefilter.")
      confirmed_paths_by_org[[org]] <- character(0)
      next
    }

    message("Organism: ", org, " - confirming compounds in KGML")
    org_hits <- character(0)

    for (path in org_paths) {
      message("  Checking ", path)

      # Build KGML URL for this pathway
      file <- paste0("https://rest.kegg.jp/get/", path, "/kgml")

      # Download KGML content
      pathway_txt <- tryCatch(
        RCurl::getURL(file),
        error = function(e) NULL
      )

      if (is.null(pathway_txt) || nchar(pathway_txt) == 0) {
        message("    ", path, " - KGML download failed, skipping")
        next
      }

      # Pause between requests
      Sys.sleep(0.5)

      # Parse KGML into KEGGPathway object
      path_parsed <- tryCatch(
        KEGGgraph::parseKGML(pathway_txt),
        error = function(e) NULL
      )

      if (is.null(path_parsed)) {
        message("    ", path, " - parseKGML error, skipping")
        next
      }

      # Convert to graph representation
      path_network <- tryCatch(
        KEGGgraph::KEGGpathway2Graph(
          path_parsed,
          genesOnly   = FALSE,
          expandGenes = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(path_network)) {
        message("    ", path, " - KEGGpathway2Graph error, skipping")
        next
      }

      # Ensure the graph contains edges
      urlcheck <- tryCatch(
        graph::edges(path_network),
        error = function(e) NULL
      )

      if (is.null(urlcheck)) {
        message("    ", path, " - no edges / bad graph, skipping")
        next
      }

      # Extract node IDs and check compound presence
      node_ids <- graph::nodes(path_network)
      if (any(all_comp_variants %in% node_ids)) {
        org_hits <- c(org_hits, path)
      }
    }

    confirmed_paths_by_org[[org]] <- unique(org_hits)
  }

  ## Return result
  confirmed_paths_by_org
}

#' Build compound-centred ego subnetworks across KEGG organisms
#'
#' Builds a compound-centred ego subnetwork around one or more KEGG
#' compounds across one or several organisms. The resulting network is
#' returned as a MetaPathNet-style edge list (3-column matrix) that can be
#' used directly for visualisation or further network analysis.
#'
#' @param compounds Character vector of KEGG compound identifiers. Only IDs in
#'   the formats \code{"Cxxxxx"} or \code{"cpd:Cxxxxx"} are accepted; any other
#'   format results in an error.
#' @param organism_codes Character vector of one or more KEGG organism codes
#'   (e.g. \code{"hsa"}, \code{"eco"}, \code{"bsu"}). For each organism, the
#'   function identifies pathways containing the input compounds and builds an
#'   ego subnetwork restricted to those pathways.
#' @param path_type Character; type of pathways to include when building the
#'   organism-specific networks. One of:
#'   \itemize{
#'     \item \code{"metabolic"}: use only metabolic pathways.
#'     \item \code{"signaling"}: use only signaling pathways.
#'     \item \code{"integrated"}: use both metabolic and signaling pathways.
#'   }
#'   The classification is derived from \code{\link{MPN_getPathIDs}}.
#' @param order Integer (default: \code{2}). Maximum graph distance (number of
#'   edges) used in \code{\link{MPN_egoNetwork}} to define the ego neighbourhood
#'   around the seed compounds. Must be a single positive integer (\eqn{\ge} 1).
#' @param mode Character; directionality of the ego search, passed to
#'   \code{\link{MPN_egoNetwork}} and ultimately \code{igraph::ego}. One of:
#'   \itemize{
#'     \item \code{"out"} (default): nodes reachable downstream from the seeds.
#'     \item \code{"in"}: nodes from which the seeds are reachable.
#'     \item \code{"all"}: treat edges as undirected for distance calculations.
#'   }
#'
#' @details
#' Organism-specific ego subnetworks are merged, deduplicated, and returned as
#' a single combined network.
#'
#' Internally, compound–pathway mapping is performed via KEGG REST / KGML
#' queries (using a non-exported helper), then routed through
#' \code{\link{MPN_keggNetwork}} and \code{\link{MPN_egoNetwork}}. These details
#' are handled automatically and do not require user intervention.
#'
#' @return
#' A character matrix with three columns:
#' \describe{
#'   \item{source}{Source node of each edge in the combined ego subnetwork.}
#'   \item{target}{Target node of each edge in the combined ego subnetwork.}
#'   \item{interaction_type}{Interaction or edge type, as returned by
#'         \code{\link{MPN_keggNetwork}} (e.g. \code{"k_compound:reversible"},
#'         \code{"k_activation"}, etc.).}
#' }
#' If no ego subnetworks can be generated for the given inputs, the function
#' stops with an error.
#'
#' @examples
#' \dontrun{
#' ## Requires extensive KEGG API querying; runtime depends on server response
#' ## 1) Define Trp/kynurenine/serotonin/anthranilate seeds
#' seed_compounds <- c(
#'   "cpd:C00078",  # L-tryptophan
#'   "cpd:C00328",  # L-kynurenine
#'   "cpd:C00780",  # serotonin
#'   "cpd:C00108"   # anthranilate
#' )
#'
#' ## 2) Build ego network (order = 2, metabolic pathways, human only)
#' ego_trp_hsa <- MPN_compoundEgoNetwork(
#'   compounds      = seed_compounds,
#'   organism_codes = "hsa",
#'   path_type      = "metabolic",
#'   order          = 2,
#'   mode           = "out"
#' )
#'
#' ## 3) (Optional) Visualise in Cytoscape
#' ## Cytoscape must be open and running
#' MPN_viewNetworkCy(
#'   network_table    = ego_trp_hsa,
#'   name             = TRUE,
#'   category         = TRUE,
#'   network_title    = "Human metabolic ego network",
#'   collection_title = "MetaPathNet_Examples"
#' )
#' }
#'
#' @seealso
#' \code{\link{MPN_getPathIDs}},
#' \code{\link{MPN_keggNetwork}},
#' \code{\link{MPN_egoNetwork}}
#'
#' @export
MPN_compoundEgoNetwork <- function(
    compounds,
    organism_codes,
    path_type = c("metabolic", "signaling", "integrated"),
    order     = 2,
    mode      = c("out", "all", "in")
) {
  #==============================================================
  # Step 0 - Basic input checks
  #==============================================================
  # Validate arguments and enforce KEGG compound ID formats
  if (missing(compounds) || length(compounds) == 0L) {
    stop("You must supply at least one compound ID in 'compounds'.")
  }
  if (missing(organism_codes) || length(organism_codes) == 0L) {
    stop("You must supply at least one KEGG organism code in 'organism_codes'.")
  }

  compounds <- unique(trimws(as.character(compounds)))
  compounds <- compounds[!is.na(compounds) & compounds != ""]

  organism_codes <- unique(trimws(as.character(organism_codes)))
  organism_codes <- organism_codes[!is.na(organism_codes) & organism_codes != ""]

  # Only accept KEGG compound IDs: Cxxxxx or cpd:Cxxxxx
  valid_comp <- grepl("^(C\\d{5}|cpd:C\\d{5})$", compounds)
  if (!all(valid_comp)) {
    bad <- paste(compounds[!valid_comp], collapse = ", ")
    stop("Invalid compound ID format detected: ", bad,
         ". Only 'Cxxxxx' or 'cpd:Cxxxxx' formats are accepted.")
  }

  path_type <- match.arg(path_type)
  mode      <- match.arg(mode)

  # Validate order for ego network
  if (length(order) != 1L || is.na(order) || !is.numeric(order) ||
      order < 1 || order != as.integer(order)) {
    stop("'order' must be a single positive integer (>= 1).")
  }
  order <- as.integer(order)

  #==============================================================
  # Step 1 - Map compounds to organism-specific pathways
  #==============================================================
  # Identify, for each organism, the KEGG pathways containing the input compounds.
  pathway_list <- MPN_mapCompoundsToPathways(
    compounds      = compounds,
    organism_codes = organism_codes
  )

  #==============================================================
  # Step 2 - Retrieve pathway annotation (metabolic / signaling)
  #==============================================================
  # Obtain pathway type information for the relevant organisms.
  path_info <- suppressMessages(
    MPN_getPathIDs(organism_codes)
  )

  #==============================================================
  # Step 3 - Standardise seed compounds to cpd: form
  #==============================================================
  # Convert compounds to "cpd:Cxxxxx" format for use as ego nodes.
  compounds_std <- unique(as.character(compounds))
  core_ids <- toupper(gsub("^(cpd:)?C", "", compounds_std, ignore.case = TRUE))
  comp_C   <- paste0("C", core_ids)
  seed_nodes <- paste0("cpd:", comp_C)

  #==============================================================
  # Step 4 - Build networks per organism and extract ego subnetworks
  #==============================================================
  # For each organism, build a pathway-specific network and extract the
  # ego neighbourhood around the seed compounds.
  ego_results <- list()

  for (org in organism_codes) {
    org_paths <- pathway_list[[org]]

    if (is.null(org_paths) || length(org_paths) == 0L) {
      message("Organism: ", org,
              " has no confirmed pathways containing the input compounds.")
      next
    }

    # Subset pathway annotation to this organism and these IDs
    org_info <- subset(
      path_info,
      Species == org & Path_id %in% org_paths
    )

    if (nrow(org_info) == 0L) {
      message("Organism: ", org,
              " has no annotated metabolic/signaling pathways for these compounds.")
      next
    }

    # Separate pathway IDs by type
    metabo_paths    <- unique(org_info$Path_id[org_info$Path_type == "metabolic"])
    signaling_paths <- unique(org_info$Path_id[org_info$Path_type == "signaling"])

    # Filter based on user choice
    if (path_type == "metabolic") {
      selected_metabo <- metabo_paths
      selected_signal <- NULL
    } else if (path_type == "signaling") {
      selected_metabo <- NULL
      selected_signal <- signaling_paths
    } else if (path_type == "integrated") {
      selected_metabo <- metabo_paths
      selected_signal <- signaling_paths
    }

    total_selected <- sum(length(selected_metabo), length(selected_signal))
    if (total_selected == 0L) {
      message("Organism: ", org, " has no ", path_type, " pathways after filtering.")
      next
    }

    message("Organism: ", org, " - ", total_selected,
            " pathways selected (",
            length(selected_metabo), " metabolic, ",
            length(selected_signal), " signaling).")

    # Build organism-specific network
    message("Building network for organism: ", org)
    net_org <- suppressMessages(
      MPN_keggNetwork(
        metabo_paths    = selected_metabo,
        signaling_paths = selected_signal
      )
    )

    # Extract ego neighbourhood around seed compounds
    message("Extracting ego neighbourhood for organism: ", org)
    ego_org <- suppressMessages(
      MPN_egoNetwork(
        network_table = net_org,
        nodes         = seed_nodes,
        order         = order,
        mode          = mode,
        output        = "matrix"
      )
    )

    ego_results[[org]] <- ego_org
  }

  #==============================================================
  # Step 5 - Combine organism-specific ego subnetworks
  #==============================================================
  # Merge per-organism ego subnetworks into a single MetaPathNet-style matrix.
  ego_results <- Filter(Negate(is.null), ego_results)

  if (length(ego_results) == 0L) {
    stop("No ego subnetworks could be generated for the given inputs.")
  }

  combined_df <- do.call(rbind, ego_results)
  combined_df <- unique(combined_df)
  rownames(combined_df) <- NULL

  # Ensure MetaPathNet-style 3-column matrix output
  combined_df <- combined_df[, c("source", "target", "interaction_type")]
  combined_mat <- as.matrix(combined_df)

  return(combined_mat)
}
