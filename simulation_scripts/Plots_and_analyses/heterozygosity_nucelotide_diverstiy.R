library(dplyr)
library(tidyr)
library(stringr)

base_dir <- "/mnt/project/CWD_reindeer/maholmen_phd/simulations/SCENARIOS/test7"
scenario_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)


# Create empty lists
all_Ho <- list()
all_pi <- list()


for (scenario_dir in scenario_dirs) {

  scenario <- basename(scenario_dir)
  message("Processing ", scenario)

  rds_file <- list.files(
    scenario_dir,
    pattern = "^Simulation_results_.*\\.rds$",
    full.names = TRUE
  )

  if (length(rds_file) == 0) next

  Simulation_results <- readRDS(rds_file[1])

  # Extract
  Ho_results <- Simulation_results$Ho_results
  pi_results <- Simulation_results$pi_results


  # Add scenario column ✅
  Ho_results$Scenario <- scenario
  pi_results$Scenario <- scenario

  # Store in list
  all_Ho[[scenario]] <- Ho_results
  all_pi[[scenario]] <- pi_results
}

# Combine everything into two datasets ✅
Ho_combined <- bind_rows(all_Ho)
pi_combined <- bind_rows(all_pi)

df <- Ho_combined
prefix <- "Ho"


##################### function ###############
process_metric <- function(df, prefix, output_name) {

  # -------------------------------
  # 1. Long format
  # -------------------------------
  long_df <- df %>%
    pivot_longer(
      cols = -c(rep, stage, Scenario),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(!is.na(value))

  # -------------------------------
  # 2. Rename scenarios
  # -------------------------------
  long_df$Scenario <- recode(long_df$Scenario,
                              "S1"  = "None",
                              "S2"  = "WT Hom ♂",
                              "S3"  = "WT All ♂",
                              "S4"  = "WT Hom",
                              "S5"  = "WT All",
                              "S6"  = "↓♂",
                              "S7"  = "WT → recover",
                              "S8"  = "↓N",
                              "S13" = "↓N = 116",
                              "S14" = "WT+225Y Hom ♂"
  )

  scenario_order <- c("None", "WT Hom ♂", "WT All ♂", "WT Hom", "WT All","↓♂", "↓N", "↓N = 116", "WT → recover", "WT+225Y Hom ♂")

  long_df$Scenario <- factor(long_df$Scenario, levels = scenario_order)

  # -------------------------------
  # 3. Filter unwanted metrics
  # -------------------------------
  #long_df <- long_df %>%
  #  filter(!grepl("control", metric, ignore.case = TRUE)) %>%
  #  filter(!grepl("SEL", metric, ignore.case = TRUE))

  # -------------------------------
  # 4. Clean metric names
  # -------------------------------
  long_df <- long_df %>%
    mutate(
      metric = case_when(
        metric == paste0(prefix, "_PRNP") ~ "ORF",
        str_starts(metric, paste0(prefix, "_PRNP_")) ~
          str_remove(metric, paste0("^", prefix, "_PRNP_")),
        str_starts(metric, paste0(prefix, "_")) ~
          str_remove(metric, paste0("^", prefix, "_"))),
      metric = ifelse(metric == "locus", "Locus", metric)
    )

  metric_order <- c("ORF", "Locus", "50Kb", "100Kb", "500Kb", "1Mb", "2Mb",
                    "CDS2", "RASSF2", "ZMYND11", "PCNA", "RPLP1", "TMEM230",
                    "PROKR2", "LARP4B", "SRPK1", "GTPBP4", "PRND", "GPCPD1",
                    "ADARB2", "SLC23A2", "DIP2C", "WDR37", "Control_25Kb",
                    "Control_500Kb", "Control_2Mb", "SEL")
  long_df$metric <- factor(long_df$metric, levels = metric_order)

  # Ensure correct stage order
  long_df$stage <- factor(long_df$stage, levels = c("start", "end"))

  # -------------------------------
  # 5. Wide format + change
  # -------------------------------
  wide_df <- long_df %>%
    group_by(rep, Scenario, metric, stage) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = value) %>%
    mutate(change = end - start)

  # -------------------------------
  # 6. T-tests
  # -------------------------------
  t_test_results <- wide_df %>%
    group_by(Scenario, metric) %>%
    summarise(
      p_value = if (n() > 1) t.test(start, end, paired = TRUE)$p.value else NA_real_,
      .groups = "drop"
    )

  # -------------------------------
  # 7. Change summary
  # -------------------------------
  change_summary <- wide_df %>%
    group_by(Scenario, metric) %>%
    summarise(
      mean_change = mean(change, na.rm = TRUE),
      sd_change   = sd(change, na.rm = TRUE),
      .groups = "drop"
    )

  # -------------------------------
  # 8. Combine
  # -------------------------------
  final_table <- change_summary %>%
    left_join(t_test_results, by = c("Scenario", "metric")) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )

  # -------------------------------
  # 9. Scale (×10⁻³)
  # -------------------------------
  final_table <- final_table %>%
    mutate(
      mean_change = mean_change * 1000,
      sd_change   = sd_change * 1000
    )

  # -------------------------------
  # 10. Format
  # -------------------------------
  table_formatted <- final_table %>%
    mutate(
      display = paste0(
        sprintf("%.2f", mean_change), " ± ",
        sprintf("%.2f", sd_change), " ",
        significance
      )
    )

  # -------------------------------
  # 11. Pivot
  # -------------------------------
  table_pivot <- table_formatted %>%
    select(Scenario, metric, display) %>%
    pivot_wider(names_from = metric, values_from = display)

  # -------------------------------
  # 12. Save
  # -------------------------------
  write.csv(
    table_pivot,
    file = file.path(base_dir, output_name),
    row.names = FALSE
  )

  return(table_pivot)
}

# For Ho
table_Ho <- process_metric(Ho_combined, prefix = "Ho",
                           output_name = "final_summary_table_HO.csv")

# For pi
table_pi <- process_metric(pi_combined, prefix = "pi",
                           output_name = "final_summary_table_PI.csv")
