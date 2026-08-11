library(RandomFieldsUtils)
library(HaploBlocker)
library(RColorBrewer)
library(vcfR)

vcf <- read.vcfR("/mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/BEAGLE_PHASED/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.CHR11_phased.beagle.out.vcf.gz")
bp_map <- getPOS(vcf)

#blocklist <- block_calculation("/mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/BEAGLE_PHASED/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.CHR11_phased.beagle.out.vcf.gz", bp_map = bp_map)
blocklist <- block_calculation("/mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/BEAGLE_PHASED/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.CHR11_phased.beagle.out.vcf.gz", bp_map = bp_map, node_min = 3, edge_min = 3)
saveRDS(blocklist, "blocklist.rds")

png("haploblocker_plot.png", width=1200, height=800, res=150)
plot_block(blocklist)
dev.off()

### find blocks around PNRP ####
my_position <- 44486053 # In the middle of the PRNP gene

find_blocks_full <- function(blocklist, position){
  results <- data.frame(
    block = integer(),
    start_bp = numeric(),
    end_bp = numeric(),
    n_haplotypes = integer(),
    haplotypes = character()
  )

  for(i in 1:length(blocklist)){
    start <- blocklist[[i]][[2]]$bp
    end <- blocklist[[i]][[3]]$bp

    if(position >= start & position <= end){
      results <- rbind(results, data.frame(
        block = i,
        start_bp = start,
        end_bp = end,
        n_haplotypes = blocklist[[i]][[5]],
        haplotypes = paste(names(blocklist[[i]][[6]]), collapse=";")
      ))
    }
  }
  return(results)
}

# Run and save
result <- find_blocks_full(blocklist, my_position)
write.csv(result, "blocks_at_position.csv", row.names=FALSE)

# find which samples this is ##
sample_names <- colnames(extract.gt(vcf))

n_samples <- length(sample_names)

haplo_map <- data.frame(
  haplo = paste0("haplo", 1:(n_samples*2)),
  sample = rep(sample_names, each=2),
  chromosome = rep(c(1,2), n_samples)
)
write.csv(haplo_map, "haplo_sample_mapping.csv", row.names=FALSE)


my_block_numbers <- result$block
#Subset the blocklist
sub_blocklist <- blocklist[my_block_numbers]

saveRDS(sub_blocklist, "sub_blocklist.rds")

# Get offset values
snp_offset <- min(sapply(sub_blocklist, function(x) x[[2]]$snp)) - 1
window_offset <- min(sapply(sub_blocklist, function(x) x[[2]]$window)) - 1

cat("SNP offset:", snp_offset, "\n")
cat("Window offset:", window_offset, "\n")

# Modify each block
sub_blocklist_mod <- sub_blocklist
for(i in seq_along(sub_blocklist_mod)){
  # Reset window numbers
  sub_blocklist_mod[[i]][[2]]$window <- sub_blocklist_mod[[i]][[2]]$window - window_offset
  sub_blocklist_mod[[i]][[3]]$window <- sub_blocklist_mod[[i]][[3]]$window - window_offset

  # Reset SNP numbers
  sub_blocklist_mod[[i]][[2]]$snp <- sub_blocklist_mod[[i]][[2]]$snp - snp_offset
  sub_blocklist_mod[[i]][[3]]$snp <- sub_blocklist_mod[[i]][[3]]$snp - snp_offset

  # Trim element 4 to new window range
  max_window <- max(sapply(sub_blocklist_mod, function(x) x[[3]]$window))
  sub_blocklist_mod[[i]][[4]] <- sub_blocklist_mod[[i]][[4]][1:max_window]
}

# Save the modified subsetted blocklist
saveRDS(sub_blocklist_mod, "sub_blocklist_mod.rds")


# Plot
png("PRNP_blocks_plot.png", width=1200, height=800, res=150)
plot_block(sub_blocklist_mod, type="snp")
dev.off()


# Load with correct column names
positions <- read.table("whole_chr11_positions.txt", header=FALSE,
                        col.names=c("snp_index", "bp_position"))

# Verify it loaded correctly
head(positions)
dim(positions)

# Get the offset that was applied
snp_offset <- min(sapply(sub_blocklist, function(x) x[[2]]$snp)) - 1

# Get SNP range from modified blocklist
min_snp_mod <- min(sapply(sub_blocklist_mod, function(x) x[[2]]$snp))
max_snp_mod <- max(sapply(sub_blocklist_mod, function(x) x[[3]]$snp))

cat("Modified SNP range:", min_snp_mod, "-", max_snp_mod, "\n")
cat("SNP offset:", snp_offset, "\n")

# Filter positions
sub_positions <- positions[positions$snp_index >= (min_snp_mod + snp_offset) &
                           positions$snp_index <= (max_snp_mod + snp_offset), ]

# Reset SNP index
sub_positions$snp_modified <- sub_positions$snp_index - snp_offset

# Reorder columns
sub_positions <- sub_positions[, c("snp_modified", "snp_index", "bp_position")]
colnames(sub_positions) <- c("snp_modified", "snp_original", "bp_position")

# Save
write.csv(sub_positions, "snp_bp_mapping_modified.csv", row.names=FALSE)
head(sub_positions)
