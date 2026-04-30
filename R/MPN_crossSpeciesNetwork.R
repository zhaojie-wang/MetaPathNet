#' Build Cross-Species MetaPathNet Network
#'
#' @description
#' Constructs and merges *MetaPathNet tables* (integrated metabolic and/or signaling networks)
#' across one or multiple KEGG organisms.
#'
#' @param organism_codes Character vector of KEGG organism codes
#'   (e.g., "hsa" for *Homo sapiens*, "eco" for *E. coli*, "bsu" for *Bacillus subtilis*).
#' @param path_type Character. One of:
#'   \describe{
#'     \item{"metabolic"}{Build network from metabolic pathways only.}
#'     \item{"signaling"}{Build network from signaling pathways only.}
#'     \item{"integrated"}{Combine both metabolic and signaling pathways into a single network.}
#'   }
#'
#' @return
#' A data frame representing the merged *MetaPathNet network*, with three columns:
#' \describe{
#'   \item{source}{Source node (KO, compound, or gene).}
#'   \item{target}{Target node (KO, compound, or gene).}
#'   \item{interaction_type}{Type of biological interaction or reaction.}
#' }
#'
#' @details
#' - KEGG pathway lists are retrieved via \code{MPN_getPathIDs()} for all requested organisms.
#' - Each species is processed separately, generating a species-specific MetaPathNet table
#'   via \code{MPN_keggNetwork()} according to \code{path_type}.
#' - Species-level tables are then merged and deduplicated into a single cross-species network.
#'
#' @examples
#' ## Cross-species metabolic network (human + E. coli strains)
#' ## Note: Network construction can be time-consuming and depends on KEGG server response
#'
#' host_microbe_code <- c("ece", "eco", "ecoo", "ecc", "ecs", "hsa")
#'
#' \dontrun{
#' ## Requires extensive KEGG API querying; runtime depends on server response
#' net_hsa_ecoli <- MPN_crossSpeciesNetwork(
#'   organism_codes = host_microbe_code,
#'   path_type      = "metabolic"   # metabolic layer only
#' )
#' }
#'
#' @seealso
#' \code{\link{MPN_getPathIDs}} for retrieving pathway lists,
#' and \code{\link{MPN_keggNetwork}} for building individual pathway networks.
#'
#' @export
MPN_crossSpeciesNetwork <- function(
    organism_codes,
    path_type = c("metabolic", "signaling", "integrated")
) {
  #==============================================================
  # Step 1 - Validate inputs and retrieve KEGG pathway metadata
  #==============================================================
  path_type <- match.arg(path_type)
  stopifnot(is.character(organism_codes), length(organism_codes) >= 1)

  message(" Retrieving KEGG pathway information across species...")
  all_paths_df <- suppressMessages(suppressWarnings(
    MPN_getPathIDs(organism_codes)
  ))

  ## Report organisms with no retrieved pathways
  requested_species <- unique(as.character(organism_codes))
  retrieved_species <- unique(all_paths_df$Species)
  missing_species   <- setdiff(requested_species, retrieved_species)

  if (length(missing_species) > 0) {
    for (sp in missing_species) {
      message(sprintf("[%s] No valid KEGG pathways retrieved. Skipping.", sp))
    }
  }

  #==============================================================
  # Step 2 - Loop over species and select pathways by type
  #==============================================================
  species_list <- unique(all_paths_df$Species)
  all_matrices <- list()

  for (sp in species_list) {
    message(sprintf("\n[%s] Processing organism...", sp))
    sp_paths <- subset(all_paths_df, Species == sp)

    # Separate pathway IDs by type
    metabo_ids    <- unique(sp_paths$Path_id[sp_paths$Path_type == "metabolic"])
    signaling_ids <- unique(sp_paths$Path_id[sp_paths$Path_type == "signaling"])

    # Filter based on user choice
    if (path_type == "metabolic") {
      selected_metabo <- metabo_ids
      selected_signal <- NULL
    } else if (path_type == "signaling") {
      selected_metabo <- NULL
      selected_signal <- signaling_ids
    } else if (path_type == "integrated") {
      selected_metabo <- metabo_ids
      selected_signal <- signaling_ids
    }

    #============================================================
    # Step 3 - Build species-specific MetaPathNet table
    #============================================================
    total_selected <- sum(length(selected_metabo), length(selected_signal))
    if (total_selected == 0) {
      message(sprintf("[%s] No %s pathways found. Skipping.", sp, path_type))
      next
    }

    message(sprintf("[%s] %d pathways selected (%d metabolic, %d signaling).",
                    sp, total_selected, length(selected_metabo), length(selected_signal)))

    message(sprintf("[%s] Building MetaPathNet table...", sp))
    mat <- tryCatch(
      suppressMessages(suppressWarnings(
        MPN_keggNetwork(
          metabo_paths    = selected_metabo,
          signaling_paths = selected_signal
        )
      )),
      error = function(e) NULL
    )

    #============================================================
    # Step 4 - Store per-species network and log status
    #============================================================
    success <- !is.null(mat)
    if (success) {
      all_matrices[[sp]] <- mat
      message(sprintf("[%s] MetaPathNet table generated successfully (%d edges).",
                      sp, nrow(mat)))
    } else {
      message(sprintf("[%s] MetaPathNet table generation failed.", sp))
    }

    ## Pause between organisms to respect KEGG REST API rate limits
    if (sp != species_list[length(species_list)]) {
      Sys.sleep(10)
    }
  }
  #==============================================================
  # Step 5 - Merge all species networks into a combined csNetwork
  #==============================================================
  if (length(all_matrices) == 0) {
    stop("No MetaPathNet tables were successfully generated for the given organisms and path_type.")
  }

  combined <- do.call(rbind, all_matrices)
  combined <- unique(combined)

  #==============================================================
  # Step 6 - Normalise rownames and report summary
  #==============================================================
  rownames(combined) <- seq_len(nrow(combined))

  message("Combined MetaPathNet network successfully built across all species.")
  message(sprintf("Total edges: %d | Total species: %d/%d",
                  nrow(combined), length(all_matrices), length(requested_species)))

  return(combined)
}
