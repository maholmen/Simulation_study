set_age_and_sex <- function(StartPop, n_total) {

  # Prepeare to set age and sex
  n_total      <- n_total
  n_calves     <- round(n_total*0.4187) # age 0 (1800 if start is 4300)
  n_adultF     <- round(n_total*0.4651) # age 2–7 females only (2000 if start is 4300)
  n_yearlings  <- n_total - n_calves - n_adultF  # age 1, minus to make sure that the number is correct
                                                 # (500 if start is 4300)

  # set age
  age_vec <- c(
    rep(0, n_calves),
    rep(1, n_yearlings),
    sample(2:10, n_adultF, replace = TRUE)
  )

  age_vec <- sample(age_vec)
  StartPop@misc$age <- age_vec

  # set sex accoring to production scheme (assuming a summer population)
  sex_vec <- character(nInd(StartPop))

  # Calves (age == 0): exactly 50/50
  is_calf <- StartPop@misc$age == 0
  n_calf  <- sum(is_calf)

  nF_calf <- floor(n_calf / 2)
  nM_calf <- n_calf - nF_calf

  sex_vec[is_calf] <- sample(
    c(rep("F", nF_calf), rep("M", nM_calf))
  )

  # Yearlings (age == 1): exactly 50/50
  is_yearling <- StartPop@misc$age == 1
  n_yearling  <- sum(is_yearling)

  nF_year <- floor(n_yearling / 2)
  nM_year <- n_yearling - nF_year

  sex_vec[is_yearling] <- sample(
    c(rep("F", nF_year), rep("M", nM_year))
  )

  # Adults (age >= 2): females only
  sex_vec[StartPop@misc$age >= 2] <- "F"

  # assign sex to population
  StartPop@sex <- sex_vec

  return(StartPop)
}

saveRDS(set_age_and_sex, file = "set_age_and_sex.rds")
