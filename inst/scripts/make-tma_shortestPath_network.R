## ==============================================================
## Script: make-tma_shortestPath_network.R
## Purpose: document how inst/extdata/tma_shortestPath_network.png
##          was generated for the MetaPathNet vignette
## ==============================================================

## Notes
## -----
## This script documents the generation of the origin-annotated
## shortest-path network figure used in the TMA/TMAO vignette.
##
## Output file:
##   inst/extdata/tma_shortestPath_network.png
##
## The shortest-path subnetwork is generated in R with MetaPathNet,
## exported to Cytoscape for visualization, and the final PNG is
## manually refined and saved to the path above.

## ==============================================================
## Step 1 — Load required packages
## ==============================================================
library(MetaPathNet)

## ==============================================================
## Step 2 — Construct the host-microbe metabolic network
## ==============================================================
network_hostMicrobe <- MPN_crossSpeciesNetwork(
  organism_codes = c("hsa", "eci"),
  path_type      = "metabolic"
)

## ==============================================================
## Step 3 — Add the manually curated FMO-associated reaction
## ==============================================================
fmo_reaction <- data.frame(
  reaction_id = "TMA_to_TMAO_FMO",
  substrates  = "cpd:C00565",
  products    = "cpd:C01104",
  ko          = "K00485",
  direction   = "reversible",
  stringsAsFactors = FALSE
)

fmo_extension <- MPN_customReaction(
  reaction_table = fmo_reaction,
  substrate_col  = "substrates",
  product_col    = "products",
  ko_col         = "ko",
  direction_col  = "direction"
)

## ==============================================================
## Step 4 — Add MetaCyc-derived reaction extensions
## ==============================================================
tma_extension <- MPN_mapReaction(
  reaction_ids = c("RXN-12900", "RXN-13946"),
  source       = "metacyc"
)

network_mixed <- MPN_mergeNetworks(
  network_hostMicrobe,
  tma_extension,
  fmo_extension
)

## ==============================================================
## Step 5 — Define selected compounds and KOs
## ==============================================================
compound_tma <- c("cpd:C00114", "cpd:C00487", "cpd:C00719",
                  "cpd:C01181", "cpd:C00565", "cpd:C01104")

gene_tma <- c("K18277", "K07811", "K00108", "K14156", "K00485", "K20038")

compound_tma <- MPN_findMappedNodes(
  nodes         = compound_tma,
  network_table = network_mixed
)$mapped_nodes

gene_tma <- MPN_findMappedNodes(
  nodes         = gene_tma,
  network_table = network_mixed
)$mapped_nodes

## ==============================================================
## Step 6 — Extract the shortest-path subnetwork
## ==============================================================
network_shortestPath <- MPN_shortestPaths(
  network_table      = network_mixed,
  source_nodes       = gene_tma,
  target_nodes       = compound_tma,
  mode               = "all",
  output             = "network_matrix",
  name               = FALSE,
  distance_threshold = 12,
  betweenness        = TRUE
)

## ==============================================================
## Step 7 — Export the origin-annotated network to Cytoscape
## ==============================================================
## Cytoscape must be open and running.

MPN_annotateOrigin(
  network_table     = network_shortestPath,
  bacteria_codes    = "eci",
  name              = TRUE,
  export_cytoscape  = TRUE,
  network_title     = "Origin-annotated shortest-path network",
  collection_title  = "MetaPathNet_Examples"
)

## ==============================================================
## Step 8 — Manual figure export
## ==============================================================
## In Cytoscape:
##   1. Open the exported network
##   2. Apply the layout used for the vignette figure
##   3. Adjust node and edge display if needed
##   4. Export the image as PNG
##
## Save the final figure to:
##   inst/extdata/tma_shortestPath_network.png

## ==============================================================
## Step 9 — Additional manual figure refinement
## ==============================================================
## The final PNG used in the vignette was manually refined after
## Cytoscape export.
##
## Manual refinements included:
##   - adding the category labels:
##       "Compounds"
##       "Microbial functions"
##       "Host functions"
##       "Shared functions"
##   - adjusting colors and layout for readability
##   - final polishing of the figure for vignette display
##
## The figure is stored in inst/extdata so that the vignette does not
## depend on Cytoscape image generation at build time. This keeps the
## vignette output stable and improves readability.
