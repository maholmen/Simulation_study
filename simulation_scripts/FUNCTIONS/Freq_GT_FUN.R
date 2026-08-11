Freq_GT <- function(
    pop,
    gen,
    rep,
    PRNP_GT,
    PRNP_Allele,
    SP
) {


  # ---------------------------------------------------------
  # PRNP genotype (diplotype) frequencies
  # ---------------------------------------------------------

  PRNP_genotypes <- pop@misc$PRNP_GT

  genotype_counts <- table(PRNP_genotypes)
  genotype_freq   <- prop.table(genotype_counts)

  PRNP_GT <- rbind(
    PRNP_GT,
    data.frame(
      rep        = rep,
      gen        = gen,
      genotype   = names(genotype_counts),
      count      = as.integer(genotype_counts),
      freq       = as.numeric(genotype_freq),
      stringsAsFactors = FALSE
    )
  )

  # ---------------------------------------------------------
  # PRNP allele frequencies
  # ---------------------------------------------------------

  alleles_from_diplo <- c(pop@misc$PRNP_Hap1, pop@misc$PRNP_Hap2)
  allele_counts <- table(alleles_from_diplo)
  allele_freq   <- prop.table(allele_counts)

  PRNP_Allele <- rbind(
    PRNP_Allele,
    data.frame(
      rep        = rep,
      gen        = gen,
      allele     = names(allele_counts),
      count      = as.integer(allele_counts),
      freq       = as.numeric(allele_freq),
      stringsAsFactors = FALSE
    )
  )


  list(
    PRNP_GT      = PRNP_GT,
    PRNP_Allele  = PRNP_Allele)
}


saveRDS(Freq_GT, file = "Freq_GT.rds")
