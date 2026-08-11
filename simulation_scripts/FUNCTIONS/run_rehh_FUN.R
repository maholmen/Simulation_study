run_rehh <- function(pop, SP, replicate_id) {

  # step 1 - extract haplotypes
  haps <- pullSegSiteHaplo(pop, simParam = SP)

  # step 2 - get SNP positions
  pos <- as.numeric(colnames(haps))

  colnames(haps) <- paste0("SNP_", seq_len(ncol(haps)))

  # Create unique haplotype IDs
  hap_ids <- paste0("hap_", rownames(haps))
  haps_with_id <- cbind(
    hap_id = hap_ids,
    haps
  )

  # convert to dataframe
  map_df <- data.frame(
    SNP = paste0("SNP_", seq_along(pos)),
    CHR = 11,
    POS = pos
  )

  write.table(
    haps_with_id,
    file = "hap_file.hap",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  write.table(
    map_df,
    file = "map.inp",
    row.names = FALSE,
    col.names = FALSE
  )

  haplohh <- data2haplohh(
    hap_file = "hap_file.hap",
    map_file = "map.inp",
    allele_coding = "01",
    chr.name = 11
  )

  # find correct markers
  marker_4 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 44485596]
  marker_385 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 44485977]
  marker_505 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 44486097]
  marker_526 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 44486118]
  marker_674 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 44486266]

  control_1_1 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 20095970]
  control_1_2 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 20109518]
  control_1_3 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 20031832]
  control_1_4 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 20033539]

  control_2_1 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 70109660]
  control_2_2 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 70182807]
  control_2_3 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 70131832]
  control_2_4 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 70012296]

  # 4 highest markers in selection area
  sel_1 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 28485308]
  sel_2 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 30237810]
  sel_3 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 28093628]
  sel_4 <- map_df$SNP[map_df$CHR == 11 & map_df$POS == 72676461] # also highly selected area


  markers <- c(
    marker_4, marker_385, marker_505, marker_526, marker_674,
    control_1_1, control_1_2, control_1_3, control_1_4,
    control_2_1, control_2_2, control_2_3, control_2_4,
    sel_1, sel_2, sel_3, sel_4
  )

  markers_id <- c(
    "marker_4", "marker_385", "marker_505", "marker_526", "marker_674",
    "control_1_1", "control_1_2", "control_1_3", "control_1_4",
    "control_2_1", "control_2_2", "control_2_3", "control_2_4",
    "sel_1", "sel_2", "sel_3", "sel_4"
  )


  #ehh_list <- lapply(markers, function(mrk) {
   # scan <- calc_ehh(haplohh, mrk = mrk)

  # this is new: remove it not working
    ehh_list <- lapply(markers, function(mrk) {

      scan <- calc_ehh(haplohh, mrk = mrk)

    data.frame(
      replicate     = replicate_id,
      marker        = mrk,
      ihh_ancestral = scan$ihh["IHH_A"],
      ihh_derived   = scan$ihh["IHH_D"]
    )
  })

  ehh_df <- do.call(rbind, ehh_list)
  ehh_df$ID <- markers_id

  return(ehh_df)
}

saveRDS(run_rehh, file = "run_rehh.rds")
