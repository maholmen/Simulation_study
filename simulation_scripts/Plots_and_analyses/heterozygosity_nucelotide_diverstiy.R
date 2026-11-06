# This code was developed with substantial assistance from OpenAI's ChatGPT. The code was reviewed and adapted by the author.

library(dplyr)
library(tidyr)
library(stringr)

base_dir <- "/mnt/project/CWD_reindeer/maholmen_phd/simulations/SCENARIOS/test7"

scenario_dirs <- list.dirs(base_dir,
                           recursive = FALSE,
                           full.names = TRUE)

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
  
  Ho_results <- Simulation_results$Ho_results
  pi_results <- Simulation_results$pi_results
  
  Ho_results$Scenario <- scenario
  pi_results$Scenario <- scenario
  
  all_Ho[[scenario]] <- Ho_results
  all_pi[[scenario]] <- pi_results
}

Ho_combined <- bind_rows(all_Ho)
pi_combined <- bind_rows(all_pi)

process_metric <- function(df, prefix, output_prefix) {
  
  # -------------------------------
  # Long format
  # -------------------------------
  long_df <- df %>%
    pivot_longer(
      cols = -c(rep, stage, Scenario),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  # -------------------------------
  # Rename scenarios
  # -------------------------------
  long_df$Scenario <- recode(long_df$Scenario,
                             "S1"  = "None",
                             "S2"  = "WT Hom ♂",
                             "S3"  = "WT All ♂",
                             "S4"  = "WT Hom ♂♀",
                             "S5"  = "WT All ♂♀",
                             "S6"  = "↓♂",
                             "S7"  = "WT → recover",
                             "S8"  = "↓N",
                             "S13" = "↓N = 116",
                             "S14" = "WT+225Y Hom ♂"
  )
  
  scenario_order <- c(
    "None",
    "WT Hom ♂",
    "WT All ♂",
    "WT Hom ♂♀",
    "WT All ♂♀",
    "↓♂",
    "↓N",
    "↓N = 116",
    "WT → recover",
    "WT+225Y Hom ♂"
  )
  
  long_df$Scenario <- factor(long_df$Scenario,
                             levels = scenario_order)
  
  # -------------------------------
  # Metric names
  # -------------------------------
  long_df <- long_df %>%
    mutate(
      metric = case_when(
        metric == paste0(prefix, "_PRNP") ~ "ORF",
        str_starts(metric, paste0(prefix, "_PRNP_")) ~
          str_remove(metric, paste0("^", prefix, "_PRNP_")),
        str_starts(metric, paste0(prefix, "_")) ~
          str_remove(metric, paste0("^", prefix, "_"))
      ),
      metric = ifelse(metric == "locus",
                      "Locus",
                      metric)
    )
  
  metric_order <- c(
    "ORF","Locus","50Kb","100Kb","500Kb","1Mb","2Mb",
    "CDS2","RASSF2","ZMYND11","PCNA","RPLP1","TMEM230",
    "PROKR2","LARP4B","SRPK1","GTPBP4","PRND","GPCPD1",
    "ADARB2","SLC23A2","DIP2C","WDR37",
    "Control_25Kb","Control_500Kb","Control_2Mb","SEL"
  )
  
  long_df$metric <- factor(long_df$metric,
                           levels = metric_order)
  
  long_df$stage <- factor(long_df$stage,
                          levels = c("start","end"))
  
  # -------------------------------
  # Wide format
  # -------------------------------
  wide_df <- long_df %>%
    group_by(rep, Scenario, metric, stage) %>%
    summarise(value = mean(value, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = stage,
                values_from = value) %>%
    mutate(
      change = end - start,
      rel_change = ((end - start) / start) * 100
    )
  
  # -------------------------------
  # T-tests
  # -------------------------------
  t_test_results <- wide_df %>%
    group_by(Scenario, metric) %>%
    summarise(
      p_value =
        if (n() > 1)
          t.test(start, end, paired = TRUE)$p.value
      else
        NA_real_,
      .groups = "drop"
    )
  
  # -------------------------------
  # STARTING VALUES
  # -------------------------------
  # starting values
  start_summary <- wide_df %>%
    group_by(metric) %>%
    summarise(
      mean_start = mean(start, na.rm = TRUE),
      sd_start   = sd(start, na.rm = TRUE),
      .groups = "drop"
    )
  # Scale only pi
  if(prefix == "pi") {
    start_summary <- start_summary %>%
    mutate(
      mean_start = mean_start * 1000,
      sd_start = sd_start * 1000
    )
  }
  
  # format
  start_table_formatted <- start_summary %>%
    mutate(
      display = paste0(
        sprintf("%.2f", mean_start), " ± ",
        sprintf("%.2f", sd_start)
      )
    )
  
  # pivot
  start_table_pivot <- start_table_formatted %>%
    select(metric, display) %>%
    pivot_wider(names_from = metric, values_from = display)
  
  write.csv(
    start_table_pivot,
    file = file.path(base_dir, paste0(output_prefix, "_starting_values.csv")),
    row.names = FALSE
  )
  
  # -------------------------------
  # ABSOLUTE CHANGE
  # -------------------------------
  abs_summary <- wide_df %>%
    group_by(Scenario, metric) %>%
    summarise(
      mean_change = mean(change, na.rm = TRUE),
      sd_change   = sd(change, na.rm = TRUE),
      .groups = "drop"
    )
  
  # combine
  abs_table <- abs_summary %>%
    left_join(t_test_results,
              by = c("Scenario","metric")) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ "ns"
      ))
  
  # scale 
  abs_table_scale <- abs_table %>%
    mutate(mean_change = mean_change * 1000,
           sd_change   = sd_change * 1000)
  
  # format
  abs_table_formatted <- abs_table_scale %>%
    mutate(
      display = paste0(
        sprintf("%.2f", mean_change), " ± ",
        sprintf("%.2f", sd_change), " ",
        significance
      )
    )
  
  # pivot 
  abs_table_pivot <- abs_table_formatted %>%
    select(Scenario, metric, display) %>%
    pivot_wider(names_from = metric, values_from = display)
  
  write.csv(
    abs_table_pivot,
    file = file.path(base_dir, paste0(output_prefix, "_absolute_change.csv")),
    row.names = FALSE
  )
  
  # -------------------------------
  # RELATIVE CHANGE (%)
  # -------------------------------
  # relative change summary
  rel_summary <- wide_df %>%
    group_by(Scenario, metric) %>%
    summarise(
      mean_rel_change = mean(rel_change, na.rm = TRUE),
      sd_rel_change   = sd(rel_change, na.rm = TRUE),
      .groups = "drop"
    )
  
  # combine
  rel_table <- rel_summary %>%
    left_join(
      t_test_results,
      by = c("Scenario", "metric")
    ) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  # format
  rel_table_formatted <- rel_table %>%
    mutate(
      display = paste0(
        sprintf("%.2f", mean_rel_change), " ± ",
        sprintf("%.2f", sd_rel_change), " ",
        significance
      )
    )
  
  # pivot
  rel_table_pivot <- rel_table_formatted %>%
    select(Scenario, metric, display) %>%
    pivot_wider(names_from = metric, values_from = display)
  
  write.csv(
    rel_table_pivot,
    file = file.path(base_dir, paste0(output_prefix, "_relative_change.csv")),
    row.names = FALSE
  )
  
  return(list(
    start_table_pivot,
    abs_table_pivot,
    rel_table_pivot
  ))
}

# Ho
Ho_tables <- process_metric(
  Ho_combined,
  prefix = "Ho",
  output_prefix = "HO"
)

# Pi
Pi_tables <- process_metric(
  pi_combined,
  prefix = "pi",
  output_prefix = "PI"
)


## PLOTTING 
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(patchwork)
library(cowplot)

base_dir <- "/mnt/project/CWD_reindeer/maholmen_phd/simulations/SCENARIOS/test7"

# --------------------------------------------------
# READ TABLES
# --------------------------------------------------

HO <- read.csv(
  file.path(base_dir, "HO_relative_change.csv"),
  check.names = FALSE
)

PI <- read.csv(
  file.path(base_dir, "PI_relative_change.csv"),
  check.names = FALSE
)

# --------------------------------------------------
# FUNCTION TO PREPARE DATA
# --------------------------------------------------

prepare_plot_data <- function(df, metric_name){
  
  df_long <- df %>%
    pivot_longer(
      cols = c("ORF","Locus","500Kb","2Mb"),
      names_to = "Region",
      values_to = "Value"
    ) %>%
    filter(Scenario != "WT+225Y Hom ♂")
  
  df_long <- df_long %>%
    mutate(
      Mean = as.numeric(
        str_extract(
          Value,
          "^-?[0-9.]+"
        )
      ),
      
      SD = as.numeric(
        str_extract(
          Value,
          "(?<=± )[-0-9.]+"
        )
      )
    )
  
  df_long$Region <- factor(
    df_long$Region,
    levels = c(
      "ORF",
      "Locus",
      "500Kb",
      "2Mb"
    )
  )
  
  df_long$Metric <- metric_name
  
  return(df_long)
}

# --------------------------------------------------
# CREATE DATASET
# --------------------------------------------------

HO_long <- prepare_plot_data(HO, "Observed heterozygosity (Hₒ)")
PI_long <- prepare_plot_data(PI, "Nucleotide diversity (π)")


# set colors: 
scenario_cols <- c(
  "None"          = "#000000",  # black
  "WT Hom ♂"      = "#0072B2",  # blue
  "WT All ♂"      = "#D55E00",  # vermillion
  "WT Hom ♂♀"     = "#56B4E9",  # sky blue
  "WT All ♂♀"     = "#CC79A7",  # reddish purple
  "↓♂"            = "#009E73",  # bluish green
  "↓N"            = "#E69F00",  # orange
  "↓N = 116"      = "#999999",  # grey
  "WT → recover"  = "#F0E442"   # yellow
)



scenario_order <- c(
  "None",
  "WT Hom ♂",
  "WT Hom ♂♀",
  "WT All ♂",
  "WT All ♂♀",
  "↓♂",
  "↓N",
  "↓N = 116",
  "WT → recover"
)



HO_long$Region <- factor(
  HO_long$Region,
  levels = c("ORF","Locus","500Kb","2Mb"),
  labels = c(
    "CDS\n(Hₒ=0.30)",
    "Locus\n(Hₒ=0.36)",
    "±250 Kb\n(Hₒ=0.31)",
    "±1 Mb\n(Hₒ=0.32)"
  )
)

HO_long$Scenario <- factor(
  HO_long$Scenario,
  levels = scenario_order
)

PI_long$Region <- factor(
  PI_long$Region,
  levels = c("ORF","Locus","500Kb","2Mb"),
  labels = c(
    "CDS\n(π=0.007 x 10⁻³)",
    "Locus\n(π=1.01 x 10⁻³)",
    "±250 Kb\n(π=0.95 x 10⁻³)",
    "±1 Mb\n(π=1.20 x 10⁻³)"
  )
)

PI_long$Scenario <- factor(
  PI_long$Scenario,
  levels = scenario_order
)


# -----------------------------
# Hₒ plot
# -----------------------------

p_HO <- ggplot(
  HO_long,
  aes(
    x = Region,
    y = Mean,
    colour = Scenario,
    group = Scenario
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    x = "",
    y = "Relative change in Hₒ (%) over 20 years",
    title = "(a) Observed heterozygosity (Hₒ)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(
      size = 11,
      face = "plain"
    ),
    legend.position = "none"
  ) + 
  scale_colour_manual(values = scenario_cols) + 
  scale_y_continuous(
    breaks = seq(-15, 20, by = 5)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))

p_HO

# -----------------------------
# π plot
# -----------------------------
p_PI <- ggplot(
  PI_long,
  aes(
    x = Region,
    y = Mean,
    colour = Scenario,
    group = Scenario
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    x = "",
    y = "Relative change in π (%) over 20 years",
    title = "(b) Nucleotide diversity (π)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(
      size = 11,
      face = "plain"
    ),
    legend.position = "none",
    legend.key.width = unit(1.5, "cm")
  ) + 
  scale_colour_manual(values = scenario_cols) + 
  scale_y_continuous(
    breaks = seq(-15, 20, by = 5)
  ) + 
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))

p_PI



# Extract legend
legend <- get_legend(
  p_PI +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 10, margin = margin(l = -10))
    ) + 
    guides(
      colour = guide_legend(
        override.aes = list(
          linetype = 0, 
          shape = 16, 
          size = 4
        )
      )
    )
)

# Convert legend to a plot
legend_plot <- wrap_elements(legend)

combined <- (
  p_HO /
    p_PI /
    legend_plot
) +
  plot_layout(
    heights = c(1, 1, 0.10)
  )


ggsave(
  file.path(base_dir, "HO_PI_lineplot.png"),
  combined,
  width = 18,
  height = 20,
  units = "cm",
  dpi = 600
)

ggsave(
  file.path(base_dir, "HO_PI_lineplot.pdf"),
  combined,
  device = cairo_pdf,
  width = 18,
  height = 20,
  units = "cm"
)




