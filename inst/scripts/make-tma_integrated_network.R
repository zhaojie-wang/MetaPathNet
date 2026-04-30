## ==============================================================
## Script: make-tma_integrated_network.R
## Purpose: document how inst/extdata/tma_integrated_network.png
##          was generated for the MetaPathNet vignette
## ==============================================================
##
## Notes
## -----
## This script documents the generation of the integrated host–microbe
## network figure used in the TMA/TMAO vignette.
##
## Output file:
##   inst/extdata/tma_integrated_network.png
##
## The network itself is generated in R with MetaPathNet, then exported
## and visualized in Cytoscape for figure production. The final PNG is
## saved manually from Cytoscape to the path above.

## ==============================================================
## Step 1 — Load required packages
## ==============================================================
library(MetaPathNet)

## ==============================================================
## Step 2 — Construct the host–microbe metabolic network
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

network_hostMicrobe <- MPN_mergeNetworks(
  network_hostMicrobe,
  fmo_extension
)

## ==============================================================
## Step 4 — Export the network to Cytoscape
## ==============================================================
## Cytoscape (version 3.9.0 or higher) must be open and running.
## The same color settings as in the vignette were used.

MPN_viewNetworkCy(
  network_table     = network_hostMicrobe,
  category          = TRUE,
  style_interaction = TRUE,
  node_compound     = "#9DC7DD",
  node_gene         = "#9ED17B",
  network_title     = "Integrated host-microbe network",
  collection_title  = "MetaPathNet_Examples"
)

## ==============================================================
## Step 5 — Manual figure export
## ==============================================================
## In Cytoscape:
##   1. Open the exported network "Integrated host-microbe network"
##   2. Apply the layout used for the vignette figure
##   3. Adjust node/edge display if needed
##   4. Export the image as PNG
##
## Save the final figure to:
##   inst/extdata/tma_integrated_network.png
##
## The exported PNG is the file used in the vignette chunk:
##   knitr::include_graphics(
##     system.file("extdata", "tma_integrated_network.png",
##                 package = "MetaPathNet")
##   )

## ==============================================================
## Step 6 — Additional manual figure refinement
## ==============================================================
## The final PNG used in the vignette was manually refined after
## Cytoscape export.
##
## Manual refinements included:
##   - adding the category labels "Compound" and "KEGG Orthology (KO)"
##   - final layout polishing for readability in the vignette
##
## The figure is stored in inst/extdata so that the vignette does not
## depend on Cytoscape image generation at build time. This keeps the
## vignette output stable and improves readability.
