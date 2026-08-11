library(AlphaSimR)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(ggpubr)

# ---- Directories ----
base_dir <- "/mnt/project/CWD_reindeer/maholmen_phd/simulations/SCENARIOS/test7"
scenario_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

# ---- Read and combine all scenarios ----
all_data <- list()

for (scenario_dir in scenario_dirs) {

  scenario <- basename(scenario_dir)
  message("Processing ", scenario)

  rds_file <- list.files(
    scenario_dir,
    pattern = "^Simulation_results_.*\\.rds$",
    full.names = TRUE
  )

  if (length(rds_file) == 0) {
    warning("Missing file in: ", scenario_dir)
    next
  }

  Simulation_results <- readRDS(rds_file[1])

  PRNP_Allele <- Simulation_results$PRNP_Allele

  # Add scenario column
  PRNP_Allele$Scenario <- scenario

  all_data[[scenario]] <- PRNP_Allele
}

# Combine everything
combined_df <- bind_rows(all_data)

# ---- CLEAN DATA ----

# Remove unwanted scenarios
combined_df <- combined_df %>%
  filter(!Scenario %in% c("S7", "S14"))

# Remove allele X
combined_df <- combined_df %>%
  filter(allele != "X")


# Rename alleles 👇
combined_df$allele <- recode(combined_df$allele,
                             "A" = "WT",
                             "B" = "225Y",
                             "D" = "176D",
                             "E" = "2M129S169M"
)

allele_order <- c("WT", "225Y", "176D", "2M129S169M")

combined_df$allele <- factor(combined_df$allele, levels = allele_order)



# ---- RENAME SCENARIOS (PUT IT HERE ✅) ----

combined_df$Scenario <- recode(combined_df$Scenario,
                           "S1"  = "None",
                           "S2"  = "WT Hom ♂",
                           "S3"  = "WT All ♂",
                           "S4"  = "WT Hom ♀♂",
                           "S5"  = "WT All ♀♂",
                           "S6"  = "↓♂",
                           "S8"  = "↓N",
                           "S13" = "↓N = 116",
                           "S14" = "WT+225Y Hom ♂"
)

scenario_order <- c("None", "WT Hom ♂", "WT All ♂", "WT Hom ♀♂", "WT All ♀♂","↓♂", "↓N", "↓N = 116", "WT+225Y Hom ♂")


combined_df$Scenario <- factor(combined_df$Scenario, levels = scenario_order)


# =========================================================
# ✅ PLOT 1: Spaghetti + mean (facet = scenario)
# =========================================================

plot_spaghetti <- ggplot(
  combined_df,
  aes(
    x = gen,
    y = freq,
    colour = allele,
    group = interaction(allele, rep, Scenario)
  )
) +
  geom_line(alpha = 0.1) +   # faint spaghetti
  # Mean lines
  stat_summary(
    aes(group = allele),
    fun = mean,
    geom = "line",
    size = 0.8
  ) +
  facet_wrap(~ Scenario, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(
    x = "Years",
    y = "Frequency",
    colour = "",
    title = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 13),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing = unit(0.7, "lines"),
    panel.border = element_blank()
  ) +
    guides(colour = guide_legend(nrow = 1))


plot_spaghetti

# Save
ggsave(
  file.path(base_dir, "Combined_spaghetti_mean.png"),
  plot_spaghetti,
  width = 12,
  height = 14,
  dpi = 300
)

### tables
wt20 <- combined_df %>%
  filter(allele == "WT") %>%
  group_by(Scenario, rep) %>%
  summarise(
    gen_WT20 = if(any(freq < 0.20))
      min(gen[freq < 0.20])
    else
      NA_real_,
    .groups = "drop"
  )

wt20_summary <- wt20 %>%
  group_by(Scenario) %>%
  summarise(
    mean_year = mean(gen_WT20, na.rm = TRUE),
    sd_year   = sd(gen_WT20, na.rm = TRUE),
    proportion =
      mean(!is.na(gen_WT20))
  )

wt10 <- combined_df %>%
  filter(allele == "WT") %>%
  group_by(Scenario, rep) %>%
  summarise(
    gen_WT10 = if(any(freq < 0.10))
      min(gen[freq < 0.10])
    else
      NA_real_,
    .groups = "drop"
  )

wt10_summary <- wt10 %>%
  group_by(Scenario) %>%
  summarise(
    mean_year = mean(gen_WT10, na.rm = TRUE),
    sd_year   = sd(gen_WT10, na.rm = TRUE),
    proportion =
      mean(!is.na(gen_WT10))
  )

wt_change <- combined_df %>%
  filter(allele == "WT") %>%
  group_by(Scenario, rep) %>%
  summarise(
    start = freq[gen == min(gen)],
    end   = freq[gen == max(gen)],
    delta = end - start,
    .groups = "drop"
  )

wt_change_summary <- wt_change %>%
  group_by(Scenario) %>%
  summarise(
    mean_start = mean(start),
    mean_end   = mean(end),
    mean_delta = mean(delta),
    sd_delta   = sd(delta)
  )

library(purrr)

annual_decline <- combined_df %>%
  filter(allele == "WT") %>%
  group_by(Scenario, rep) %>%
  group_modify(~{

    fit <- lm(freq ~ gen, data = .x)

    tibble(
      slope = coef(fit)[2]
    )
  })

annual_decline_summary <- annual_decline %>%
  group_by(Scenario) %>%
  summarise(
    mean_decline = mean(slope),
    sd_decline   = sd(slope)
  )
