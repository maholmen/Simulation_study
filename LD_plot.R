library(ggplot2)
library(dplyr)
library(patchwork)

ld <- read.table("LD_PRNP_R2_0.1.ld", header = TRUE, stringsAsFactors = FALSE)
ld$distance <- abs(ld$BP_A - ld$BP_B)

snp_info <- data.frame(
  SNP_A = c(
    "chromosome_11_OX459947.2_557_44485977_G_A",
    "chromosome_11_OX459947.2_557_44486118_A_G",
    "chromosome_11_OX459947.2_557_44486266_A_C"
  ),
  BP = c(44485977, 44486118, 44486266),
  label = c("Codon 2, Codon 129 and Codon 169 represented by Codon 129", "Codon 176", "Codon 225")
)

make_ld_plot <- function(snp_id, bp, label, xlab_text = "") {
  
  df <- ld %>% dplyr::filter(SNP_A == snp_id)
  
  ggplot(df, aes(BP_B, R2)) +
    geom_point() +
    geom_point(aes(x = bp, y = 1.0), color = "red", size = 2) +
    xlim(43487546, 45486053) +
    labs(title = label, x = xlab_text, y = NULL) +
    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.0, size = 10),
      axis.text = element_text(size = 10)
    ) + 
    labs(y = expression(R^2))
  
}

plots <- mapply(
  make_ld_plot,
  snp_info$SNP_A,
  snp_info$BP,
  snp_info$label,
  SIMPLIFY = FALSE
)







LD_plot <- wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave(
  "LD_plot_3_panel.png",
  LD_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)


