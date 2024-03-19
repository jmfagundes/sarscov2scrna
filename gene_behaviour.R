source("source_me.R")

# read/create rank matrices

load(paste0("RData/", out.dir, "cell_infection_pseudotime_matrices_TL.RData"))

if (use.TPT) {
  hBECs.mtx <- hBECs.mtx.TPT
  colon.mtx <- colon.mtx.TPT
  ileum.mtx <- ileum.mtx.TPT
  
  hBECs.TL <- hBECs.TL.TPT
  colon.TL <- colon.TL.TPT
  ileum.TL <- ileum.TL.TPT
  rm(hBECs.mtx.TPT, colon.mtx.TPT, ileum.mtx.TPT,
     hBECs.TL.TPT, colon.TL.TPT, ileum.TL.TPT)
}

hBECs.rank <- process.rank(matrices = hBECs.mtx$matrices[c("infected", "uninfected_0dpi")], method = "random")
colon.rank <- process.rank(matrices = colon.mtx$matrices[c("infected", "uninfected_0hpi")], method = "random")
ileum.rank <- process.rank(matrices = ileum.mtx$matrices[c("infected", "uninfected_0hpi")], method = "random")

# calculate RSI and 0 rate for infected cells

hBECs.RSI <- calc.RSI(hBECs.rank["infected"])
colon.RSI <- calc.RSI(colon.rank["infected"])
ileum.RSI <- calc.RSI(ileum.rank["infected"])

hBECs.0rate <- calc.0_rate(hBECs.mtx$matrices[c("infected", "uninfected_0dpi")])
colon.0rate <- calc.0_rate(colon.mtx$matrices[c("infected", "uninfected_0hpi")])
ileum.0rate <- calc.0_rate(ileum.mtx$matrices[c("infected", "uninfected_0hpi")])

# shuffle RSI, higher dropout due to higher viral load is confounding this analysis

hBECs.RSI.boot <- calc.RSI(hBECs.rank["infected"], shuffle.rep = 1000, genes.to.keep = names(hBECs.0rate$infected[hBECs.0rate$infected < .1]))
colon.RSI.boot <- calc.RSI(colon.rank["infected"], shuffle.rep = 1000, genes.to.keep = names(colon.0rate$infected[colon.0rate$infected < .1]))
ileum.RSI.boot <- calc.RSI(ileum.rank["infected"], shuffle.rep = 1000, genes.to.keep = names(ileum.0rate$infected[ileum.0rate$infected < .1]))

# reorder RSI by accumulated expression

hBECs.RSI$infected <- hBECs.RSI$infected[match(names(sort(hBECs.mtx$matrices$infected %>% rowSums(), decreasing = TRUE)),
                                               names(hBECs.RSI$infected))]
colon.RSI$infected <- colon.RSI$infected[match(names(sort(colon.mtx$matrices$infected %>% rowSums(), decreasing = TRUE)),
                                               names(colon.RSI$infected))]
ileum.RSI$infected <- ileum.RSI$infected[match(names(sort(ileum.mtx$matrices$infected %>% rowSums(), decreasing = TRUE)),
                                               names(ileum.RSI$infected))]

# plot fit to power law

gg.fit <- plot.power_law(list(hBECs = hBECs.TL$log_log.mean_sd$infected,
                              colon = colon.TL$log_log.mean_sd$infected,
                              ileum = ileum.TL$log_log.mean_sd$infected),
                         list(hBECs = hBECs.RSI$infected,
                              colon = colon.RSI$infected,
                              ileum = ileum.RSI$infected), "RSI", c(0, 100)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = .25) +
  geom_abline(intercept = -1.7, slope = .5, color = "black", linewidth = .25, linetype = "dashed") +
  scale_color_viridis() +
  xlab("mean (log)") + ylab("standard deviation (log)")

ggsave("figs_out/final/RSI_regression.pdf", gg.fit, width = 6.85, height = 3)
ggsave(paste0("figs_out/", out.dir, "RSI_regression_simplified.png"), gg.fit, width = 6.85, height = 3, dpi = 300)
