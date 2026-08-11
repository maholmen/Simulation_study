###

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
reps <- as.numeric(args[1]) # how many replicates to run
S_haplo <- as.numeric(args[2]) # how many haplotypes to sample (might change for smaller test - for anything larger than 1000, sample 1000)
gen_expand <- as.numeric(args[3]) # over how many generations the population shall expand (most 4)
n_startPop <- as.numeric(args[4]) # how big the startingpop should be #4300 normally
gens <- as.numeric(args[5]) # how many generations to run
strategy_file <- args[6] # what to remove and how many breeding males

# set some numbers based on n_startPop
n_total      <- n_startPop
n_calves     <- round(n_total*0.4187) # age 0 (1800 if start is 4300)
n_adultF     <- round(n_total*0.4651) # age 2–7 females only (2000 if start is 4300)
n_yearlings  <- n_total - n_calves - n_adultF  # age 1, minus to make sure that the number is correct
# (500 if start is 4300)


# extract clean strategy name
strategy_name <- tools::file_path_sans_ext(basename(strategy_file))

# Make table
params <- data.frame(
  reps = reps,
  S_haplo = S_haplo,
  gen_expand = gen_expand,
  n_startPop = n_startPop,
  gens = gens,
  strategy = strategy_name,
  n_total = n_total,
  n_calves = n_calves,
  n_adultF = n_adultF,
  n_yearlings = n_yearlings
)


write.csv(params,
          file = paste0("simulation_parameters_", strategy_name, ".csv"),
          row.names = FALSE)

# load libraries to use

library(AlphaSimR)
library(rehh)
library(hierfstat)

# Load haplotype data, genmap and functions
sim_data <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/alphasim_input_data.rds")
annotatePrnp <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/FUNCTIONS/annotatePrnp.rds")
set_age_and_sex <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/FUNCTIONS/set_age_and_sex.rds")
Freq_GT <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/FUNCTIONS/Freq_GT.rds")
Ho_Pi_td <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/FUNCTIONS/Ho_Pi_td.rds")
run_rehh <- readRDS("/mnt/project/CWD_reindeer/maholmen_phd/simulations/FUNCTIONS/run_rehh.rds")

## run before every simulation ##
# define results containers
PRNP_GT <- data.frame()
PRNP_Allele <- data.frame()
Ho_results <- data.frame()
all_ehh_df <- data.frame()
start_ehh_df <- data.frame()
Ho_results <- data.frame()
pi_results <- data.frame()
td_results <- data.frame()

# create founderPop - based on actual data
founderPop <- importHaplo(
  haplo  = sim_data$haplo_bin,
  genMap = sim_data$genMap
)


# run replicates - sampling of 1000 haplotypes, expanding population and selection
for (rep in 1:reps) {
  # Sample haplotypes
  S <- sampleHaplo(founderPop, S_haplo, inbred = FALSE, ploidy = 2, replace = TRUE)

  # set global simulation parameters
  SP <- SimParam$new(S)
  SP$setSexes("yes_rand")

  # create population object
  pop <- newPop(S, simParam = SP)

  # expand population size over a few generations (4 generations give ~5000 animals)
  for (gen in 1:gen_expand) {
    parents <- pop

    offspring <- randCross(
      parents,
      nCrosses = (nInd(pop)/2), #based on assumption that approximately 50% is female
      nProgeny = 1,
      simParam = SP
    )

    # combine offspring and parents
    pop <- c(parents, offspring)
  }

  StartPop <- sample(pop, n_startPop) # rename to StartPop

  # annotate PRNP
  StartPop <- annotatePrnp(
    pop     = StartPop
  )

  # set age and sex
  StartPop <- set_age_and_sex(StartPop, n_startPop)

  Start <- StartPop[StartPop@misc$age > 0]

  # check EHH before I start
  ehh_df_start <- run_rehh(
    pop = Start,
    SP = SP,
    replicate_id = rep
  )

  start_ehh_df <- rbind(start_ehh_df, ehh_df_start)

  # check Heterozygosity, nucelotide diversity and TajimaD on StartPop
  Ho_Pi_td_result_start <- Ho_Pi_td(Start, SP, rep = rep, stage = "start")

  Ho_results <- rbind(Ho_results, Ho_Pi_td_result_start$Ho)
  pi_results <- rbind(pi_results, Ho_Pi_td_result_start$pi)
  td_results <- rbind(td_results, Ho_Pi_td_result_start$td)

  # remove unecessary data
  rm(offspring, parents, pop, S)

  # define population as another name to not overwrite StartPop
  pop <- StartPop

  # define calves before starting to fit in the loop
  C1 <-  pop[pop@misc$age == 0] # last loops calves

  for (gen in 1:gens) {
    # culling based on strategy file
    source(strategy_file)

    # cull so it is only 500 calves left, 250 bucks and 250 females
    C1XM <- (m <- C1_R[isMale(C1_R)])[sample(nInd(m), round(n_yearlings/2))] #culled male calves
    C1XF <- (f <- C1_R[isFemale(C1_R)])[sample(nInd(f), round(n_yearlings/2))] #n_yearling is defined in the set_age_and_sex function

    C1X <- c(C1XM, C1XF) # last loops calves - culled this year - this will be the breeding males next loop

    # remove so it is only 1750 females between the age of 2-10:
    AF_keep <- (isFemale(pop) & pop@misc$age <= 10 & pop@misc$age >= 2) # is female and between 2 and 10 years
    AF <- pop[AF_keep]

    # cull ~250 females (keep 1750)
    n_f <- min(round(n_adultF*0.875), nInd(AF)) #keep it flexible in case there is under the set number
    AFX <- AF[sample(nInd(AF), n_f)] # culled adult females - more control

    # breeding females is this loops yearlings and adult females (~2000 animals)
    # female yearling
    YF <- pop[isFemale(pop) & pop@misc$age == 1]

    BF <- c(AFX, YF) # this loops breeding females

    # Breeding population:
    B <- c(BF, BM)

    # create mating plan
    BF_IDs <- sample(BF@id, floor(nInd(BF)*0.9), replace = FALSE) # set to 0.9 because this is the average "kalvetilgang" after loss
    BM_IDs   <- sample(BM@id, length(BF_IDs), replace = TRUE) # this can easily be adjusted so
    # some males contribute more.

    # mating happening in the fall, but calves are not "here" until the spring
    C0 <- makeCross(
      B,
      crossPlan = cbind(BF_IDs, BM_IDs),
      nProgeny  = 1,
      simParam  = SP
    )
    # add random mutations to calves
    AlphaSimR::mutate(C0, mutRate = 1e-08, returnPos = FALSE, simParam = SP)

    # set all metadata (misc) for calves
    C0 <- annotatePrnp(pop = C0) # PRNP genotype and "A-dosage"
    C0@misc$age <- rep(0, nInd(C0)) # age

    # prepare data for next loop, the remaining population is culled calves and
    # breeding females. C0 is part of next years population.
    # combine calves and breeding females
    pop <- mergePops(list(BF, C1X))

    # age population
    pop@misc$age <- pop@misc$age + 1

    # set C0 as C1
    C1 <- C0

    # collect data
    out <- Freq_GT(
      pop = pop,
      gen = gen,
      rep = rep,
      PRNP_GT = PRNP_GT,
      PRNP_Allele = PRNP_Allele,
      SP = SP
    )

    PRNP_GT <- out$PRNP_GT
    PRNP_Allele <- out$PRNP_Allele

    # remove unnecessary objects
    rm(out, AF, AFX, B, BF, BM, C0, C1X, C1XF, C1XM, f,
       m, YF, AF_keep, BF_IDs, BM_IDs, n_f)

  }

  ehh_df <- run_rehh(
    pop = pop,
    SP = SP,
    replicate_id = rep
  )

  all_ehh_df <- rbind(all_ehh_df, ehh_df)

  Ho_Pi_td_result_end <- Ho_Pi_td(pop, SP, rep = rep, stage = "end")

  Ho_results <- rbind(Ho_results, Ho_Pi_td_result_end$Ho)
  pi_results <- rbind(pi_results, Ho_Pi_td_result_end$pi)
  td_results <- rbind(td_results, Ho_Pi_td_result_end$td)


}

Simulation_results <- list(
  PRNP_Allele = PRNP_Allele,
  PRNP_GT = PRNP_GT,
  Ho_results = Ho_results,
  all_ehh_df = all_ehh_df,
  start_ehh_df = start_ehh_df,
  pi_results = pi_results,
  td_results = td_results
)

# save with strategy name
outfile <- paste0("Simulation_results_", strategy_name, ".rds")

saveRDS(Simulation_results, outfile)
