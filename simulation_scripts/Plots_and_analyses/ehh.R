library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ---------------------------------------------------------
# DIRECTORIES
# ---------------------------------------------------------
base_dir <- "/mnt/project/CWD_reindeer/maholmen_phd/simulations/SCENARIOS/test7"
scenario_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

# ---------------------------------------------------------
# STORE ALL DATA
# ---------------------------------------------------------
all_ehh <- list()

# ---------------------------------------------------------
# LOOP THROUGH SCENARIOS
# ---------------------------------------------------------
all_ehh <- list()

for (scenario_dir in scenario_dirs) {

  scenario <- basename(scenario_dir)
  message("Processing ", scenario)

  # -------------------------
  # Load file
  # -------------------------
  rds_file <- list.files(
    scenario_dir,
    pattern = "^Simulation_results_.*\\.rds$",
    full.names = TRUE
  )

  if (length(rds_file) == 0) next

  Simulation_results <- readRDS(rds_file[1])

  # -------------------------
  # Extract EHH data
  # -------------------------
  start_ehh_df <- Simulation_results$start_ehh_df
  end_ehh_df   <- Simulation_results$all_ehh_df

  # ✅ skip if missing
  if (is.null(start_ehh_df) || is.null(end_ehh_df)) {
    message("Skipping ", scenario, ": missing EHH data")
    next
  }

  start_ehh_df$time <- "start"
  end_ehh_df$time   <- "end"

  ehh_all <- dplyr::bind_rows(start_ehh_df, end_ehh_df)

  # ✅ check required columns BEFORE pivot
  if (!all(c("ihh_ancestral", "ihh_derived") %in% colnames(ehh_all))) {
    message("Skipping ", scenario, ": missing ihh columns")
    next
  }

  # -------------------------
  # FIX derived / ancestral swap
  # -------------------------
  is_674 <- ehh_all$ID == "marker_674"
  tmp <- ehh_all$ihh_ancestral[is_674]
  ehh_all$ihh_ancestral[is_674] <- ehh_all$ihh_derived[is_674]
  ehh_all$ihh_derived[is_674]   <- tmp

  # -------------------------
  # Rename IDs
  # -------------------------
  id_map <- c(
    "control_1_1" = "C1", "control_1_2" = "C2", "control_1_3" = "C3",
    "control_1_4" = "C4", "control_1_5" = "C5",
    "control_2_1" = "C6", "control_2_2" = "C7", "control_2_3" = "C8",
    "marker_4"    = "PRNP codon 2",
    "marker_385"  = "PRNP codon 129",
    "marker_505"  = "PRNP codon 169",
    "marker_526"  = "PRNP codon 176",
    "marker_674"  = "PRNP codon 225",
    "sel_1" = "S1", "sel_2" = "S2", "sel_3" = "S3", "sel_4" = "S4"
  )

  ehh_all$ID <- id_map[ehh_all$ID]

  # -------------------------
  # Convert to long format
  # -------------------------
  ehh_long <- ehh_all %>%
    tidyr::pivot_longer(
      cols = c(ihh_ancestral, ihh_derived),
      names_to  = "allele",
      values_to = "ihh"
    )

  # Add scenario
  ehh_long$Scenario <- scenario

  all_ehh[[scenario]] <- ehh_long
}

# ---------------------------------------------------------
# COMBINE ALL SCENARIOS ✅
# ---------------------------------------------------------
ehh_long <- bind_rows(all_ehh)

ehh_reduced <- ehh_long %>%
  filter(!Scenario %in% c("S1","S2", "S3", "S4", "S5", "S6", "S7", "S8", "S14"))

ehh_long <- ehh_long %>%
  filter(!Scenario %in% c("S13", "S14"))


# ---------------------------------------------------------
# CLEAN & FORMAT
# ---------------------------------------------------------

# rename allele
ehh_long$allele <- recode(ehh_long$allele,
                          "ihh_ancestral" = "Wildtype",
                          "ihh_derived"   = "Alternative"
)

# remove NA
ehh_long <- ehh_long %>% filter(!is.na(ihh))

# ---------------------------------------------------------
# RENAME SCENARIOS (same as before)
# ---------------------------------------------------------

ehh_long$Scenario <- recode(ehh_long$Scenario,
                            "S1"  = "None",
                            "S2"  = "WT Hom ♂",
                            "S3"  = "WT All ♂",
                            "S4"  = "WT Hom ♀♂",
                            "S5"  = "WT All ♀♂",
                            "S6"  = "↓♂",
                            "S8"  = "↓N",
                            "S14" = "WT+225Y Hom ♂",
                            "S7"  = "WT → recover"
)



scenario_order <- c("None", "WT Hom ♂", "WT All ♂", "WT Hom ♀♂", "WT All ♀♂","↓♂", "↓N", "↓N = 116", "WT+225Y Hom ♂", "WT → recover")

ehh_long$Scenario <- factor(ehh_long$Scenario, levels = scenario_order)

allele_order <- c("Wildtype", "Alternative")

ehh_long$allele <- factor(ehh_long$allele, levels = allele_order)

# ---------------------------------------------------------
# CALCULATE CHANGE ✅
# ---------------------------------------------------------
ehh_wide <- ehh_long %>%
  group_by(replicate, Scenario, ID, allele, time) %>%
  summarise(ihh = mean(ihh), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = ihh) %>%
  mutate(change = end - start)


ehh_wide <- ehh_wide %>%
  mutate(
    change_rel = (end - start) / start
  )

library(dplyr)

change_tests <- ehh_wide %>%
  group_by(Scenario, ID, allele) %>%
  summarise(
    p_value = if (n() > 1) t.test(change)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )


ehh_plot_data <- ehh_wide %>%
  filter(grepl("^PRNP", ID))

ehh_plot_data <- ehh_wide %>%
  filter(grepl("^PRNP", ID)) %>%
  left_join(change_tests, by = c("Scenario", "ID", "allele"))

library(stringr)

ehh_plot_data$ID <- str_remove(ehh_plot_data$ID, "^PRNP ")
ehh_plot_data$ID <- str_replace(ehh_plot_data$ID, "^codon", "Codon")

ID_order <- c("Codon 2", "Codon 129", "Codon 169", "Codon 176", "Codon 225")
ehh_plot_data$ID <- factor(ehh_plot_data$ID, levels = ID_order)

ehh_plot_data <- ehh_plot_data %>%
  mutate(sig_group = ifelse(p_value < 0.05, "Significant", "Not significant"))


ehh_plot_data$log <- ehh_plot_data$change*0.001

EHH <- ggplot(ehh_plot_data,
       aes(x = Scenario, y = log, fill)) +
 coord_cartesian(ylim = c(-100, 160)) +
  geom_violin(trim = TRUE, alpha = 0.8, fill = "grey") +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    color = "black",
    width = 0.5,
    linewidth = 0.3
  ) +
  facet_grid(ID ~ allele) +
  geom_text(
    aes(label = significance, y = 140),
    stat = "summary",
    fun = max,
    size = 5
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    strip.text = element_text(size = 14),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12)

  ) +
  labs(
    x = "",
    y = "Δ iHH × 10⁻³",
    title = ""
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = seq(-80, 120, by = 20)
  )

EHH
ggsave(
  file.path(base_dir, "EHH.png"),
  EHH,
  width = 12,
  height = 14,
  dpi = 300
)


summary_tbl <- ehh_plot_data %>%
  group_by(Scenario, ID, allele) %>%
  summarise(
    mean_change = mean(change, na.rm = TRUE),
    sd_change = sd(change, na.rm = TRUE),
    mean_rel = mean(change_rel, na.rm = TRUE) * 100,
    p_value = first(p_value),
    significance = first(significance),
    .groups = "drop"
  )

summary_tbl <- summary_tbl %>%
  mutate(
    Delta_iHH =
      sprintf("%.1f ± %.1f",
              mean_change,
              sd_change),
    Relative_change =
      sprintf("%.1f%%",
              mean_rel)
  )

########### bottleneck scenario ###################

# ---------------------------------------------------------
# CLEAN & FORMAT
# ---------------------------------------------------------

# rename allele
ehh_reduced$allele <- recode(ehh_reduced$allele,
                          "ihh_ancestral" = "Wildtype",
                          "ihh_derived"   = "Alternative"
)

# remove NA
ehh_reduced <- ehh_reduced %>% filter(!is.na(ihh))

# ---------------------------------------------------------
# RENAME SCENARIO (same as before)
# ---------------------------------------------------------

ehh_reduced$Scenario <- recode(ehh_reduced$Scenario,
                            "S13"  = "↓N = 116")


allele_order <- c("Wildtype", "Alternative")

ehh_reduced$allele <- factor(ehh_reduced$allele, levels = allele_order)

# ---------------------------------------------------------
# CALCULATE CHANGE ✅
# ---------------------------------------------------------
ehh_wide <- ehh_reduced %>%
  group_by(replicate, Scenario, ID, allele, time) %>%
  summarise(ihh = mean(ihh), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = ihh) %>%
  mutate(change = end - start)


change_tests <- ehh_wide %>%
  group_by(Scenario, ID, allele) %>%
  summarise(
    p_value = if (n() > 1) t.test(change)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )


ehh_plot_data <- ehh_wide %>%
  filter(grepl("^PRNP", ID))

ehh_plot_data <- ehh_wide %>%
  filter(grepl("^PRNP", ID)) %>%
  left_join(change_tests, by = c("Scenario", "ID", "allele"))


ehh_plot_data$ID <- str_remove(ehh_plot_data$ID, "^PRNP ")
ehh_plot_data$ID <- str_replace(ehh_plot_data$ID, "^codon", "Codon")

ID_order <- c("Codon 2", "Codon 129", "Codon 169", "Codon 176", "Codon 225")
ehh_plot_data$ID <- factor(ehh_plot_data$ID, levels = ID_order)

ehh_plot_data <- ehh_plot_data %>%
  mutate(sig_group = ifelse(p_value < 0.05, "Significant", "Not significant"))



summary_tbl_scenario13 <- ehh_plot_data %>%
group_by(Scenario, ID, allele) %>%
summarise(
  mean_change = mean(change, na.rm = TRUE),
  sd_change = sd(change, na.rm = TRUE),
  mean_rel = mean(change_rel, na.rm = TRUE) * 100,
  p_value = first(p_value),
  significance = first(significance),
  .groups = "drop"
)

summary_tbl_scenario13 <- summary_tbl_scenario13 %>%
  mutate(
    Delta_iHH =
      sprintf("%.1f ± %.1f",
              mean_change,
              sd_change),
    Relative_change =
      sprintf("%.1f%%",
              mean_rel)
  )
