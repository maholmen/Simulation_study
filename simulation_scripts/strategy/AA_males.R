# remove all A calves - Scenario 5
remove <- (isMale(C1) & C1@misc$A_dosage == 2)
C1_R <- C1[!remove]

# breeding males is this years male yearling, but will assign them as BM
BM_1 <- pop[isMale(pop) & pop@misc$age == 1]

BM <- BM_1
