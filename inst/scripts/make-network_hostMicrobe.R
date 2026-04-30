## ==============================================================
## Script: make-network_hostMicrobe.R
## Purpose: document how inst/extdata/network_hostMicrobe.rds
##          was generated for the MetaPathNet vignette
## ==============================================================
##
## Object:
##   network_hostMicrobe.rds
##
## Biological context:
##   Host–microbiome metabolic network (human + E. coli O18:K1:H7 UTI89)
##   used in the choline–TMA–TMAO vignette.
##
## Source function:
##   MPN_crossSpeciesNetwork(organism_codes = c("hsa", "eci"),
##                           path_type = "metabolic")
##
## Note:
##   This object is precomputed because live KEGG network construction
##   can take roughly 5 minutes and depends on KEGG server availability.
##   To regenerate, source this script with internet access and
##   the MetaPathNet package installed.
## ==============================================================

library(MetaPathNet)

## Build the host–microbiome metabolic network
network_hostMicrobe <- MPN_crossSpeciesNetwork(
  organism_codes = c("hsa", "eci"),
  path_type      = "metabolic"
)

## Save to inst/extdata
out_dir <- file.path("inst", "extdata")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(network_hostMicrobe, file = file.path(out_dir, "network_hostMicrobe.rds"))

message("Done. Saved to inst/extdata/network_hostMicrobe.rds")
