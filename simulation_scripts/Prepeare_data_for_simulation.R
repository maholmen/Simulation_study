# SECTION 1: Load and prepare haplotype data
haplo <- as.matrix(read.table(
  "haplo_matrix.txt",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
))

# Sort markers by physical position
numeric_colnames <- as.numeric(colnames(haplo))
haplo_order      <- haplo[, order(numeric_colnames)]


# Force to be binary, this means that multi-alleic sites now will act as bi-allelic sites
haplo_bin <- haplo_order
haplo_bin[haplo_bin > 1] <- 1  # convert anything >1 to 1


# SECTION 2: Build genetic map (1 cM / Mb)

marker_bp  <- as.numeric(colnames(haplo_bin))

# cM positions
gen_pos <- marker_bp * 1e-8 # it says Morgans in the documentation, but this is
                            # something I am unsure about

genMap <- data.frame(
  markerName = colnames(haplo_bin),
  chromosome = 1,       # only chr11, treated as single chromosome
  position   = gen_pos
)

# Create a list of objects to save
sim_data <- list(
  genMap          = genMap,
  haplo_bin   = haplo_bin
)

# Save to disk
saveRDS(sim_data, file = "alphasim_input_data.rds")
