# keep all 
C1_R <- C1

# breeding males is this years male yearling, but will assign them as BM
BM_1 <- pop[isMale(pop) & pop@misc$age == 1]
BM <- BM_1
