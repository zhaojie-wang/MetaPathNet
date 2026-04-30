## Build example data for MetaPathNet documentation
##
## This script builds a precomputed host-microbe tryptophan-related
## network used in examples, tests, and vignettes.
##
## Run manually from the package root:
##
## devtools::load_all()
## source("data-raw/MetaPathNet_example_network.R")

if (!dir.exists("data")) {
  dir.create("data")
}

## Human tryptophan metabolism and related immune/signaling pathways
metabo_paths_hsa <- "hsa00380"
signaling_paths_hsa <- c(
  "hsa04060", "hsa04630", "hsa04064",
  "hsa04660", "hsa04659"
)

## Build human network
net_trp_hsa <- MPN_keggNetwork(
  metabo_paths    = metabo_paths_hsa,
  signaling_paths = signaling_paths_hsa
)

## Build E. coli network
net_trp_eco <- MPN_keggNetwork(
  metabo_paths    = c("eco00380", "eco00400"),
  signaling_paths = NULL
)

## Merge host and microbial networks
MetaPathNet_example_network <- MPN_mergeNetworks(
  net_trp_hsa,
  net_trp_eco
)

## Keep a matrix-format edge list for downstream MetaPathNet examples
MetaPathNet_example_network <- as.matrix(MetaPathNet_example_network)

## Save as package data
save(
  MetaPathNet_example_network,
  file = "data/MetaPathNet_example_network.rda",
  compress = "xz"
)
