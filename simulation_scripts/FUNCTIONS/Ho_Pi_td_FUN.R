Ho_Pi_td <- function(pop, SP, rep, stage) {

  ## Determine intervals
  S_ORF <- 44485593
  E_ORF <- 44486363
  C <- 20000000
  S1 <- 25000000
  S2 <- 32000000

  marker_ids <- as.numeric(names(SP$genMap[[1]]))

  regions <- list(

    PRNP = list(
      markers = names(SP$genMap[[1]])[marker_ids >= S_ORF & marker_ids <= E_ORF],
      length  = E_ORF - S_ORF
    ),

    PRNP_locus = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= 44467081 & marker_ids <= 44491590
      ],
      length = 44491590 - 44467081
    ),

    PRNP_50Kb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (S_ORF - 25000) & marker_ids <= (E_ORF + 25000)
      ],
      length = (E_ORF + 25000) - (S_ORF - 25000)
    ),

    PRNP_100Kb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (S_ORF - 50000) & marker_ids <= (E_ORF + 50000)
      ],
      length = (E_ORF + 50000) - (S_ORF - 50000)
    ),

    PRNP_500Kb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (S_ORF - 250000) & marker_ids <= (E_ORF + 250000)
      ],
      length = (E_ORF + 250000) - (S_ORF - 250000)
    ),

    PRNP_1Mb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (S_ORF - 500000) & marker_ids <= (E_ORF + 500000)
      ],
      length = (E_ORF + 500000) - (S_ORF - 500000)
    ),

    PRNP_2Mb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (S_ORF - 1000000) & marker_ids <= (E_ORF + 1000000)
      ],
      length = (E_ORF + 1000000) - (S_ORF - 1000000)
    ),

    Control_25Kb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (C - 12500) & marker_ids <= (C + 12500)
      ],
      length = (C + 12500) - (C - 12500)
    ),

    Control_500Kb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (C - 250000) & marker_ids <= (C + 250000)
      ],
      length = (C + 250000) - (C - 250000)
    ),

    Control_2Mb = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= (C - 1000000) & marker_ids <= (C + 1000000)
      ],
      length = (C + 1000000) - (C - 1000000)
    ),

    SEL = list(
      markers = names(SP$genMap[[1]])[
        marker_ids >= S1 & marker_ids <= S2
      ],
      length = S2 - S1
    ),

    CDS2 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44938083 & marker_ids <= 44953262],
      length  = 44953262 - 44938083
    ),

    RASSF2 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44560979 & marker_ids <= 44583593],
      length  = 44583593 - 44560979
    ),


    ZMYND11 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44254248 & marker_ids <= 44284596],
      length  = 44284596 - 44254248
    ),

    IDI1 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 43772249 & marker_ids <= 43774445],
      length  = 43774445 - 43772249
    ),

    PCNA = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44881746 & marker_ids <= 44886218],
      length  = 44886218 - 44881746
    ),

    RPLP1 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 45074165 & marker_ids <= 45074509],
      length  = 45074509 - 45074165
    ),


    TMEM230 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44873189 & marker_ids <= 44878810],
      length  = 44878810 - 44873189
    ),

    PROKR2 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 45026741 & marker_ids <= 45034118],
      length  = 45034118 - 45026741
    ),

    LARP4B = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 43835395 & marker_ids <= 43914307],
      length  = 43914307 - 43835395
    ),

    SRPK1 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44913236 & marker_ids <= 44913667],
      length  = 44913667 - 44913236
    ),

    GTPBP4 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 43778455 & marker_ids <= 43792446],
      length  = 43792446 - 43778455
    ),

    PRND = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44513238 & marker_ids <= 44513774],
      length  = 44513774 - 44513238
    ),


    GPCPD1 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 45242182 & marker_ids <= 45273292],
      length  = 45273292 - 45242182
    ),

    ADARB2 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 43665268 & marker_ids <= 43718310],
      length  = 43718310 - 43665268
    ),

    SLC23A2 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44644276 & marker_ids <= 44726965],
      length  = 44726965 - 44644276
    ),

    DIP2C = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 44190216 & marker_ids <= 44235914],
      length  = 44235914 - 44190216
    ),

    WDR37 = list(
      markers = names(SP$genMap[[1]])[marker_ids >= 43743122 & marker_ids <= 43758615],
      length  = 43758615 - 43743122
    )
  )

  # 1. extract dosages
  dosages <- lapply(regions, function(region) {
    pullMarkerGeno(
      pop,
      markers  = region$markers,
      simParam = SP
    )
  })

  # 2. Ho
  Ho_list <- lapply(dosages, function(dos) {
    mean(colMeans(dos == 1, na.rm = TRUE), na.rm = TRUE)
  })

  # 3. pi
  pi_list <- lapply(names(dosages), function(name) {
    pi.dosage(
      dosages[[name]],
      L = regions[[name]]$length
    )
  })
  names(pi_list) <- names(dosages)

  # 4. TajimaD
  td_list <- lapply(dosages, function(dos) {
    TajimaD.dosage(dos)
  })

  # 5. convert to one-row dataframes
  Ho_df <- setNames(as.data.frame(as.list(Ho_list)), paste0("Ho_", names(Ho_list)))
  pi_df <- setNames(as.data.frame(as.list(pi_list)), paste0("pi_", names(pi_list)))
  td_df <- setNames(as.data.frame(as.list(td_list)), paste0("td_", names(td_list)))

  # 6. return everything together
  list(
    Ho = cbind(rep = rep, stage = stage, Ho_df),
    pi = cbind(rep = rep, stage = stage, pi_df),
    td = cbind(rep = rep, stage = stage, td_df)
  )
}

saveRDS(Ho_Pi_td, file = "Ho_Pi_td.rds")
