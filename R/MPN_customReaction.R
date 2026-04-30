#' Build MetaPathNet-style edges from user-defined reactions
#'
#' Convert one or more manually curated reactions into the standard
#' MetaPathNet 3-column edge list format
#' (\code{source}, \code{target}, \code{interaction_type}).
#'
#' Each input row represents one reaction. Substrates and products can contain
#' one or multiple node IDs, separated by semicolons. Optional KO annotations
#' can be provided to build compound-KO-compound edges; otherwise the function
#' returns direct compound-compound edges. User-defined compounds or enzymes can
#' first be mapped to KEGG identifiers with \code{MPN_keggFinder()} when possible.
#'
#' @param reaction_table A data.frame in which each row corresponds to one
#'   custom reaction.
#' @param substrate_col Character. Name of the column containing substrate node
#'   IDs. Multiple IDs must be separated by semicolons.
#' @param product_col Character. Name of the column containing product node
#'   IDs. Multiple IDs must be separated by semicolons.
#' @param ko_col Optional character. Name of the column containing KO IDs.
#'   Multiple IDs must be separated by semicolons. If \code{NULL}, no KO layer
#'   is added and direct substrate-product edges are returned.
#' @param direction_col Optional character. Name of the column containing
#'   reaction direction. Accepted values are \code{"reversible"},
#'   \code{"irreversible"}, and \code{"unknown"}. If \code{NULL}, all reactions
#'   are treated as \code{"unknown"}.
#'
#' @return
#' A character matrix with three columns:
#' \describe{
#'   \item{source}{Source node ID.}
#'   \item{target}{Target node ID.}
#'   \item{interaction_type}{Reaction type, e.g.
#'         \code{"custom:reversible"}, \code{"custom:irreversible"}, or
#'         \code{"custom:unknown"}.}
#' }
#'
#' @examples
#' ## Example — Add a user-defined irreversible reaction
#'
#' ## 1) Build a custom reaction table
#' ##    Multiple substrates, products, or KO terms must be separated by semicolons.
#' custom_df <- data.frame(
#'   reaction_id = "custom_reaction_1",
#'   substrates  = "cpd:C00078; cpd:C00001",
#'   products    = "cpd:C00328",
#'   ko          = "K00453",
#'   direction   = "irreversible",
#'   stringsAsFactors = FALSE
#' )
#'
#' ## 2) Convert the user-defined reaction into MetaPathNet-style edges
#' custom_edges <- MPN_customReaction(
#'   reaction_table = custom_df,
#'   substrate_col  = "substrates",
#'   product_col    = "products",
#'   ko_col         = "ko",
#'   direction_col  = "direction"
#' )
#'
#' custom_edges
#'
#' @seealso
#' \code{\link{MPN_keggFinder}}, \code{\link{MPN_mapReaction}},
#' \code{\link{MPN_mergeNetworks}}
#'
#' @export
MPN_customReaction <- function(
    reaction_table,
    substrate_col = "substrates",
    product_col   = "products",
    ko_col        = NULL,
    direction_col = NULL
) {
  #==============================================================
  # Step 1 — Basic checks and input normalisation
  #==============================================================
  if (missing(reaction_table) || is.null(reaction_table)) {
    stop("reaction_table must be provided.")
  }

  rt <- as.data.frame(reaction_table, stringsAsFactors = FALSE)

  if (nrow(rt) == 0) {
    stop("reaction_table must contain at least one reaction.")
  }

  required_cols <- c(substrate_col, product_col)
  if (!all(required_cols %in% colnames(rt))) {
    stop(
      "reaction_table must contain the specified substrate/product columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  if (!is.null(ko_col) && !(ko_col %in% colnames(rt))) {
    stop("ko_col was provided but is not present in reaction_table.")
  }

  if (!is.null(direction_col) && !(direction_col %in% colnames(rt))) {
    stop("direction_col was provided but is not present in reaction_table.")
  }

  #==============================================================
  # Step 2 — Loop over reactions and build edge list
  #==============================================================
  all_edges <- list()
  edge_index <- 1L

  for (i in seq_len(nrow(rt))) {
    subs_raw <- rt[[substrate_col]][i]
    prods_raw <- rt[[product_col]][i]

    if (is.na(subs_raw) || is.na(prods_raw) || subs_raw == "" || prods_raw == "") {
      next
    }

    subs <- trimws(unlist(strsplit(as.character(subs_raw), ";", fixed = TRUE)))
    prods <- trimws(unlist(strsplit(as.character(prods_raw), ";", fixed = TRUE)))

    subs <- unique(subs[subs != ""])
    prods <- unique(prods[prods != ""])

    if (length(subs) == 0 || length(prods) == 0) {
      next
    }

    if (!is.null(direction_col)) {
      direction <- tolower(trimws(as.character(rt[[direction_col]][i])))
      if (is.na(direction) || direction == "") direction <- "unknown"
    } else {
      direction <- "unknown"
    }

    if (!(direction %in% c("reversible", "irreversible", "unknown"))) {
      warning("Skipping reaction row ", i, ": invalid direction value '", direction, "'.")
      next
    }

    itype <- paste0("custom:", direction)

    # Optional KO layer
    if (!is.null(ko_col)) {
      ko_raw <- rt[[ko_col]][i]

      if (!is.na(ko_raw) && ko_raw != "") {
        kos <- trimws(unlist(strsplit(as.character(ko_raw), ";", fixed = TRUE)))
        kos <- unique(kos[kos != ""])
      } else {
        kos <- character(0)
      }
    } else {
      kos <- character(0)
    }

    #==============================================================
    # Step 3 — Build edges for this reaction
    #==============================================================
    if (length(kos) > 0) {
      # Substrates -> KO
      for (ko in kos) {
        for (sub in subs) {
          all_edges[[edge_index]] <- data.frame(
            source = sub,
            target = ko,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source = ko,
              target = sub,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }

        # KO -> products
        for (prod in prods) {
          all_edges[[edge_index]] <- data.frame(
            source = ko,
            target = prod,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source = prod,
              target = ko,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }
      }
    } else {
      # Direct substrate -> product edges
      for (sub in subs) {
        for (prod in prods) {
          all_edges[[edge_index]] <- data.frame(
            source = sub,
            target = prod,
            interaction_type = itype,
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L

          if (direction == "reversible") {
            all_edges[[edge_index]] <- data.frame(
              source = prod,
              target = sub,
              interaction_type = itype,
              stringsAsFactors = FALSE
            )
            edge_index <- edge_index + 1L
          }
        }
      }
    }
  }

  #==============================================================
  # Step 4 — Merge all reactions and return clean edge matrix
  #==============================================================
  if (length(all_edges) == 0) {
    warning("No valid custom reactions could be converted.")
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
