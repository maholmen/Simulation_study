library(ggplot2)
library(dplyr)
library(patchwork)

# Read LD data
ld <- read.table(
  "LD_PRNP_R2_0.1.ld",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Common plotting function
make_ld_plot <- function(snp_id, bp, label) {
  
  df <- ld %>%
    filter(SNP_A == snp_id)
  
  ggplot(df, aes(BP_B, R2)) +
    geom_point(size = 0.8) +
    geom_point(
      aes(x = bp, y = 1),
      colour = "red",
      size = 2
    ) +
    xlim(43487546, 45486053) +
    labs(
      title = label,
      x = NULL,
      y = expression(R^2)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0, size = 10),
      axis.text = element_text(size = 10)
    )
}

############################################################
# 5-PANEL FIGURE
############################################################

snp_info_5 <- data.frame(
  SNP_A = c(
    "chromosome_11_OX459947.2_557_44485596_G_A",
    "chromosome_11_OX459947.2_557_44485977_G_A",
    "chromosome_11_OX459947.2_557_44486097_G_A",
    "chromosome_11_OX459947.2_557_44486118_A_G",
    "chromosome_11_OX459947.2_557_44486266_A_C"
  ),
  BP = c(
    44485596,
    44485977,
    44486097,
    44486118,
    44486266
  ),
  label = c(
    "(a) Codon 2",
    "(b) Codon 129",
    "(c) Codon 169",
    "(d) Codon 176",
    "(e) Codon 225"
  )
)

plots_5 <- mapply(
  make_ld_plot,
  snp_info_5$SNP_A,
  snp_info_5$BP,
  snp_info_5$label,
  SIMPLIFY = FALSE
)

LD_plot_5 <- wrap_plots(plots_5, ncol = 1) &
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave(
  "LD_plot_5_panel.png",
  LD_plot_5,
  width = 18,
  height = 24,
  units = "cm",
  dpi = 600
)

ggsave(
  "LD_plot_5_panel.pdf",
  LD_plot_5,
  device = cairo_pdf,
  width = 18,
  height = 24,
  units = "cm"
)

############################################################
# 3-PANEL FIGURE
############################################################

snp_info_3 <- data.frame(
  SNP_A = c(
    "chromosome_11_OX459947.2_557_44485977_G_A",
    "chromosome_11_OX459947.2_557_44486118_A_G",
    "chromosome_11_OX459947.2_557_44486266_A_C"
  ),
  BP = c(
    44485977,
    44486118,
    44486266
  ),
  label = c(
    "(a) Codons 2, 129 and 169 (represented by Codon 129)",
    "(b) Codon 176",
    "(c) Codon 225"
  )
)

plots_3 <- mapply(
  make_ld_plot,
  snp_info_3$SNP_A,
  snp_info_3$BP,
  snp_info_3$label,
  SIMPLIFY = FALSE
)

LD_plot_3 <- wrap_plots(plots_3, ncol = 1) &
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave(
  "LD_plot_3_panel.png",
  LD_plot_3,
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

ggsave(
  "LD_plot_3_panel.pdf",
  LD_plot_3,
  device = cairo_pdf,
  width = 18,
  height = 16,
  units = "cm"
)
