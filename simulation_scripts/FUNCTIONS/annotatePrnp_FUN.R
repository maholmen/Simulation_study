annotatePrnp <- function(pop) {

  profiles <- list(
    A = c(0, 0, 0, 0, 1),
    B = c(0, 0, 0, 0, 0),
    D = c(0, 0, 0, 1, 1),
    E = c(1, 1, 1, 0, 1)
  )

  # The 5 coding variant positions — stored as character to match column names
  profiles_positions <- as.character(c(44485596, 44485977, 44486097, 44486118, 44486266))

  # ── pull haplotypes from AlphaSimR ───────────────────────────────────────
  prnp_haps <- pullMarkerHaplo(pop, markers = profiles_positions, haplo = "all", asRaw = FALSE, simParam = NULL)
  # dimensions: individuals × markers × haplotypes(2)

  nInd <- dim(prnp_haps)[1]

  # ── exact haplotype assignment (your logic) ───────────────────────────────
  assign_profile <- function(x, profiles) {
    for (nm in names(profiles)) {
      if (all(x == profiles[[nm]])) return(nm)
    }
    return("X")
  }

  hap_labels <- apply(prnp_haps, 1, assign_profile, profiles = profiles)

  hap1 = hap_labels[seq(1, length(hap_labels), by = 2)]
  hap2 = hap_labels[seq(2, length(hap_labels), by = 2)]

  # ── Genotypes (A/B not B/A) ──────────────────────────────────────────────
  genotypes <- sapply(
    seq(1, length(hap_labels), by = 2),
    function(i) {
      paste(sort(c(hap_labels[i], hap_labels[i+1])), collapse = "/")
    }
  )

  # ── store everything in misc ──────────────────────────────────────────────
  pop@misc$PRNP_Hap1        <- hap1
  pop@misc$PRNP_Hap2        <- hap2
  pop@misc$PRNP_GT          <- genotypes

  # assign A_dosage
  pop@misc$A_dosage <- (pop@misc$PRNP_Hap1 == "A") + (pop@misc$PRNP_Hap2 == "A") # A-dosage
  # assign B_dosage
  pop@misc$B_dosage <- (pop@misc$PRNP_Hap1 == "B") + (pop@misc$PRNP_Hap2 == "B") # B-dosage
  # assign D_dosage
  pop@misc$D_dosage <- (pop@misc$PRNP_Hap1 == "D") + (pop@misc$PRNP_Hap2 == "D") # D-dosage
  # assign E_dosage
  pop@misc$E_dosage <- (pop@misc$PRNP_Hap1 == "E") + (pop@misc$PRNP_Hap2 == "E") # E-dosage

  return(pop)
 }

saveRDS(annotatePrnp, file = "annotatePrnp.rds")
