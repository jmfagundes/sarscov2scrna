source("source_me.R")

# get target number of cell when building the simulated dataset

load("RData/MT/cell_infection_pseudotime_matrices_TL.RData")

hBECs.sim.target.counts <- hBECs.mtx$matrices$infected %>% colSums()
colon.sim.target.counts <- colon.mtx$matrices$infected %>% colSums()
ileum.sim.target.counts <- ileum.mtx$matrices$infected %>% colSums()

# load data

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

dir.create("figs_out/final")

# bin pseudotimes

hBECs.bin30.TL <- pseudotime_bin(hBECs.mtx$matrices$infected, n.bins = hBECs.n.bins)
colon.bin30.TL <- pseudotime_bin(colon.mtx$matrices$infected, n.bins = colon.n.bins)
ileum.bin30.TL <- pseudotime_bin(ileum.mtx$matrices$infected, n.bins = ileum.n.bins)

pseudotime_bin30.TL_params <- rbind(cbind(hBECs.bin30.TL$params, cell = "hBECs"),
                                    cbind(colon.bin30.TL$params, cell = "colon"),
                                    cbind(ileum.bin30.TL$params, cell = "ileum"))

pseudotime_bin30.TL_params$cell <- factor(pseudotime_bin30.TL_params$cell, levels = c("hBECs", "colon", "ileum"))

# ANOVA

V.bin30.aov <- aov(lm(V ~ pseudotime_bin * cell + sparsity + n.genes, data = pseudotime_bin30.TL_params))
beta.bin30.aov <- aov(lm(beta ~ pseudotime_bin * cell + sparsity + n.genes, data = pseudotime_bin30.TL_params))

V.bin30.pes <- V.bin30.aov %>% effectsize::eta_squared()
beta.bin30.pes <- beta.bin30.aov %>% effectsize::eta_squared()

# plot power laws

hBECs.llmsd.bin30 <- setNames(hBECs.bin30.TL$log_log.mean_sd, nm = 1:length(hBECs.bin30.TL$matrices)) %>% bind_rows(.id = "bin")
colon.llmsd.bin30 <- setNames(colon.bin30.TL$log_log.mean_sd, nm = 1:length(colon.bin30.TL$matrices)) %>% bind_rows(.id = "bin")
ileum.llmsd.bin30 <- setNames(ileum.bin30.TL$log_log.mean_sd, nm = 1:length(ileum.bin30.TL$matrices)) %>% bind_rows(.id = "bin")

hBECs.llmsd.bin30$bin <- as.numeric(hBECs.llmsd.bin30$bin)
colon.llmsd.bin30$bin <- as.numeric(colon.llmsd.bin30$bin)
ileum.llmsd.bin30$bin <- as.numeric(ileum.llmsd.bin30$bin)

hBECs.llmsd.bin30$gene <- rownames(hBECs.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)
colon.llmsd.bin30$gene <- rownames(colon.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)
ileum.llmsd.bin30$gene <- rownames(ileum.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)

hBECs.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", hBECs.llmsd.bin30$gene), "MT", "nuclear")
colon.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", colon.llmsd.bin30$gene), "MT", "nuclear")
ileum.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", ileum.llmsd.bin30$gene), "MT", "nuclear")

# annotate some genes in hBECs that always follow an exponential distribution

# to check manually, in particular t = 8, 9, 10, 12, 14, 18, 20, 28

# ggplot(hBECs.llmsd.bin30 %>% mutate(gene = rownames(hBECs.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)), aes(x = mean, y = sd, color = is.mt, label = gene)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
#  geom_point(size = .01) +
#  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
#  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black")) + geom_text_repel(size = .5)
# ggsave("figs_out/TEMP_hBECs_30bin_ggrepel.pdf", width = 50, height = 50, limitsize = FALSE)

hBECs.nuclear.exponential.genes <- c("SCGB3A1", "SCGB1A1", "S100A2", "SERPINB3", "WFDC2", "SLPI", "S100A9", "LCN2") 

# add exponential nuclear genes for hBECs

hBECs.llmsd.bin30$is.mt[hBECs.llmsd.bin30$gene %in% hBECs.nuclear.exponential.genes] <- "exp.nuc"

# plot

hBECs.llmsd.bin30.gg <- ggplot(hBECs.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("orange", "green", "black"))
colon.llmsd.bin30.gg <- ggplot(colon.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("colon") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))
ileum.llmsd.bin30.gg <- ggplot(ileum.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("ileum") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))

# export figures

only.bins.gg <- ggplot(pseudotime_bin30.TL_params, aes(x = V, y = beta, color = pseudotime_bin)) +
  labs(color = "pseudotime") +
  geom_point(size = 1.2) + scale_shape_manual(values = c(NA, 18)) +
  facet_wrap(~cell) + geom_path(linewidth = .25) +
  guides(shape = "none")

ggsave("figs_out/TPT_MT/taylors_params_30_pseudotime_bin.pdf", only.bins.gg + 
         theme(axis.text = element_text(size = 6), legend.title = element_text(size = 8), legend.text = element_text(size = 6)), 
       width = 6.85, height = 3.5)

ggsave("figs_out/TPT_MT/hBECs_fit_30_pseudotime_bin.png", hBECs.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85, dpi = 300)
ggsave("figs_out/TPT_MT/colon_fit_30_pseudotime_bin.png", colon.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85, dpi = 300)
ggsave("figs_out/TPT_MT/ileum_fit_30_pseudotime_bin.png", ileum.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85, dpi = 300)

ggsave("figs_out/TPT_MT/hBECs_fit_30_pseudotime_bin.pdf", hBECs.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85)
ggsave("figs_out/TPT_MT/colon_fit_30_pseudotime_bin.pdf", colon.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85)
ggsave("figs_out/TPT_MT/ileum_fit_30_pseudotime_bin.pdf", ileum.llmsd.bin30.gg +
         theme(axis.text = element_text(size = 4.25), legend.title = element_text(size = 8), legend.text = element_text(size = 6)),
       width = 6.85, height = 6.85)

# save data

save(hBECs.bin30.TL,
     colon.bin30.TL,
     ileum.bin30.TL,
     file = "RData/TPT_MT/cell_infection_pseudotime_bin_matrices_TL.RData")

# fit segmented TL
# remove mitochondrial genes since they behave differently and filter genes at .7 sparsity to avoid low expressing exponential genes

hBECs.mtx.TPT.noMT <- process.gene.count("Expected_TPT_matrix_hBECs_final.csv", remove.MT = TRUE, TPT = TRUE, infected.cells = names(hBECs.mtx$viral_accumulation$infected))
colon.mtx.TPT.noMT <- process.gene.count("Expected_TPT_matrix_colon_final.csv", remove.MT = TRUE, TPT = TRUE, infected.cells = names(colon.mtx$viral_accumulation$infected))
ileum.mtx.TPT.noMT <- process.gene.count("Expected_TPT_matrix_ileum_final.csv", remove.MT = TRUE, TPT = TRUE, infected.cells = names(ileum.mtx$viral_accumulation$infected))

# remove some exponential nuclear genes from hBECs

hBECs.mtx.TPT.noMT$matrices$infected <- hBECs.mtx.TPT.noMT$matrices$infected[!rownames(hBECs.mtx.TPT.noMT$matrices$infected) %in% hBECs.nuclear.exponential.genes,]

hBECs.bin30.TL.TPT.noMT <- pseudotime_bin(hBECs.mtx.TPT.noMT$matrices$infected, zero.rate.threshold = .7, n.bins = hBECs.n.bins)
colon.bin30.TL.TPT.noMT <- pseudotime_bin(colon.mtx.TPT.noMT$matrices$infected, zero.rate.threshold = .7, n.bins = colon.n.bins)
ileum.bin30.TL.TPT.noMT <- pseudotime_bin(ileum.mtx.TPT.noMT$matrices$infected, zero.rate.threshold = .7, n.bins = ileum.n.bins)

# 1 breakpoint

hBECs.bin30.seg.TL <- fit.seg.TL(hBECs.bin30.TL.TPT.noMT)
colon.bin30.seg.TL <- fit.seg.TL(colon.bin30.TL.TPT.noMT)
ileum.bin30.seg.TL <- fit.seg.TL(ileum.bin30.TL.TPT.noMT)

breakpoint.sparsity.bin30.tb <- rbind(hBECs.bin30.seg.TL$params %>% cbind(cell = "hBECs"),
                                      colon.bin30.seg.TL$params %>% cbind(cell = "colon"),
                                      ileum.bin30.seg.TL$params %>% cbind(cell = "ileum"))

colnames(breakpoint.sparsity.bin30.tb)[1] <- "pseudotime_bin"

# check plots

hBECs.llmsd.bin30.TPT.noMT <- setNames(hBECs.bin30.TL.TPT.noMT$log_log.mean_sd, nm = 1:length(hBECs.bin30.TL.TPT.noMT$matrices)) %>% bind_rows(.id = "bin")
colon.llmsd.bin30.TPT.noMT <- setNames(colon.bin30.TL.TPT.noMT$log_log.mean_sd, nm = 1:length(colon.bin30.TL.TPT.noMT$matrices)) %>% bind_rows(.id = "bin")
ileum.llmsd.bin30.TPT.noMT <- setNames(ileum.bin30.TL.TPT.noMT$log_log.mean_sd, nm = 1:length(ileum.bin30.TL.TPT.noMT$matrices)) %>% bind_rows(.id = "bin")

hBECs.llmsd.bin30.TPT.noMT$bin <- as.numeric(hBECs.llmsd.bin30.TPT.noMT$bin)
colon.llmsd.bin30.TPT.noMT$bin <- as.numeric(colon.llmsd.bin30.TPT.noMT$bin)
ileum.llmsd.bin30.TPT.noMT$bin <- as.numeric(ileum.llmsd.bin30.TPT.noMT$bin)

hBECs.llmsd.bin30.TPT.noMT$gene <- rownames(hBECs.llmsd.bin30.TPT.noMT) %>% gsub("\\.*", "", .)
colon.llmsd.bin30.TPT.noMT$gene <- rownames(colon.llmsd.bin30.TPT.noMT) %>% gsub("\\.*", "", .)
ileum.llmsd.bin30.TPT.noMT$gene <- rownames(ileum.llmsd.bin30.TPT.noMT) %>% gsub("\\.*", "", .)

hBECs.llmsd.bin30.TPT.noMT.gg <- ggplot(hBECs.llmsd.bin30.TPT.noMT, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")
colon.llmsd.bin30.TPT.noMT.gg <- ggplot(colon.llmsd.bin30.TPT.noMT, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("colon") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")
ileum.llmsd.bin30.TPT.noMT.gg <- ggplot(ileum.llmsd.bin30.TPT.noMT, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("ileum") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")

# ANOVA and other tests

breakpoint.sparsity.bin30.tb <- rbind(hBECs.bin30.seg.TL$params %>% cbind(cell = "hBECs"),
                                      colon.bin30.seg.TL$params %>% cbind(cell = "colon"),
                                      ileum.bin30.seg.TL$params %>% cbind(cell = "ileum")) %>%
  separate(V, into = c("V1", "V2"), sep = ",") %>%
  separate(beta, into = c("beta1", "beta2"), sep = ",")

breakpoint.sparsity.bin30.tb$cell <- factor(breakpoint.sparsity.bin30.tb$cell, levels = c("hBECs", "colon", "ileum"))

colnames(breakpoint.sparsity.bin30.tb)[1] <- "pseudotime_bin"

breakpoint.sparsity.bin30.tb$pseudotime_bin <- breakpoint.sparsity.bin30.tb$pseudotime_bin %>% as.numeric()
breakpoint.sparsity.bin30.tb$breakpoint <- breakpoint.sparsity.bin30.tb$breakpoint %>% as.numeric() %>% exp() %>% log10()
breakpoint.sparsity.bin30.tb$V1 <- breakpoint.sparsity.bin30.tb$V1 %>% as.numeric()
breakpoint.sparsity.bin30.tb$V2 <- breakpoint.sparsity.bin30.tb$V2 %>% as.numeric()
breakpoint.sparsity.bin30.tb$beta1 <- breakpoint.sparsity.bin30.tb$beta1 %>% as.numeric()
breakpoint.sparsity.bin30.tb$beta2 <- breakpoint.sparsity.bin30.tb$beta2 %>% as.numeric()

# test breakpoint

breakpoint.lm <- lm(breakpoint ~ pseudotime_bin * cell + n.genes + sparsity,
                    data = breakpoint.sparsity.bin30.tb)

breakpoint.aov <- aov(breakpoint.lm)

breakpoint.pes <- effectsize::eta_squared(breakpoint.lm)

# test estimated V's and beta's with the breakpoint as a covariate

V1.lm <- lm(V1 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
            data = breakpoint.sparsity.bin30.tb)
beta1.lm <- lm(beta1 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
               data = breakpoint.sparsity.bin30.tb)

V2.lm <- lm(V2 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
            data = breakpoint.sparsity.bin30.tb)
beta2.lm <- lm(beta2 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
               data = breakpoint.sparsity.bin30.tb)

V1.aov <- aov(V1.lm)
beta1.aov <- aov(beta1.lm)
V2.aov <- aov(V2.lm)
beta2.aov <- aov(beta2.lm)

V1.pes <- effectsize::eta_squared(V1.lm)
beta1.pes <- effectsize::eta_squared(beta1.lm)
V2.pes <- effectsize::eta_squared(V2.lm)
beta2.pes <- effectsize::eta_squared(beta2.lm)

# plot

breakpoint.sparsity.bin30.gg.tb <- breakpoint.sparsity.bin30.tb %>%
  dplyr::rename(`V[1]` = V1, `V[2]` = V2, `beta[1]` = beta1, `beta[2]` = beta2, `breakpoint~(log)` = breakpoint) %>%
  dplyr::select(pseudotime_bin, `V[1]`, `V[2]`, `beta[1]`, `beta[2]`, `breakpoint~(log)`, cell) %>%
  pivot_longer(-c(pseudotime_bin, cell))

breakpoint.sparsity.bin30.gg.tb$name <- factor(breakpoint.sparsity.bin30.gg.tb$name,
                                               levels = c("breakpoint~(log)", "V[1]", "beta[1]", "V[2]", "beta[2]"))

npsi1.breakpoint.gg <- ggplot(breakpoint.sparsity.bin30.gg.tb, aes(pseudotime_bin, value, color = cell)) + geom_line() +
  facet_grid(rows = vars(name), scales = "free", switch = "y", labeller = label_parsed) + 
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + xlab("pseudotime bin")

ggsave("figs_out/final/seg_params.pdf", npsi1.breakpoint.gg, width = 6.85, height = 7)

# save data

save(hBECs.bin30.seg.TL,
     colon.bin30.seg.TL,
     ileum.bin30.seg.TL,
     file = "RData/TPT_MT/cell_infection_pseudotime_bin_matrices_seg_TL.RData")

# simulate noise in uninfected cells

# reorder vector of target counts

hBECs.sim.target.counts <- hBECs.sim.target.counts[match(names(hBECs.mtx$viral_accumulation$infected), names(hBECs.sim.target.counts))]
colon.sim.target.counts <- colon.sim.target.counts[match(names(colon.mtx$viral_accumulation$infected), names(colon.sim.target.counts))]
ileum.sim.target.counts <- ileum.sim.target.counts[match(names(ileum.mtx$viral_accumulation$infected), names(ileum.sim.target.counts))]

set.seed(1024)
hBECs.mtx.uninfected.down.sample <- down.sample(hBECs.mtx$matrices$uninfected_0dpi[sample(length(hBECs.mtx$viral_accumulation$infected))],
                                                hBECs.sim.target.counts)
colon.mtx.uninfected.down.sample <- down.sample(colon.mtx$matrices$uninfected_0hpi,
                                                colon.sim.target.counts %>% tail(ncol(colon.mtx$matrices$uninfected_0hpi)))
ileum.mtx.uninfected.down.sample <- down.sample(ileum.mtx$matrices$uninfected_0hpi,
                                                ileum.sim.target.counts %>% tail(ncol(ileum.mtx$matrices$uninfected_0hpi)))

# to simulate an infected dataset

#hBECs.mtx.uninfected.down.sample <- down.sample(hBECs.mtx$matrices$infected,
#                                                hBECs.sim.target.counts)
#colon.mtx.uninfected.down.sample <- down.sample(colon.mtx$matrices$infected,
#                                                colon.sim.target.counts)
#ileum.mtx.uninfected.down.sample <- down.sample(ileum.mtx$matrices$infected,
#                                                ileum.sim.target.counts)

hBECs.uninfected.sim.bin30.TL <- pseudotime_bin(hBECs.mtx.uninfected.down.sample, hBECs.n.bins)
colon.uninfected.sim.bin30.TL <- pseudotime_bin(colon.mtx.uninfected.down.sample, colon.n.bins)
ileum.uninfected.sim.bin30.TL <- pseudotime_bin(ileum.mtx.uninfected.down.sample, ileum.n.bins)

pseudotime_uninfected.sim.bin30.TL_params <- rbind(cbind(hBECs.uninfected.sim.bin30.TL$params, cell = "hBECs"),
                                                   cbind(colon.uninfected.sim.bin30.TL$params, cell = "colon"),
                                                   cbind(ileum.uninfected.sim.bin30.TL$params, cell = "ileum"))

pseudotime_uninfected.sim.bin30.TL_params$cell <- factor(pseudotime_uninfected.sim.bin30.TL_params$cell, levels = c("hBECs", "colon", "ileum"))

# see if zero rate per bin per gene is similar to infected cells

hBECs.bin30.0rate.tb <- calc.0_rate(hBECs.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)
colon.bin30.0rate.tb <- calc.0_rate(colon.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)
ileum.bin30.0rate.tb <- calc.0_rate(ileum.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)

hBECs.uninfected.sim.bin30.0rate.tb <- calc.0_rate(hBECs.uninfected.sim.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)
colon.uninfected.sim.bin30.0rate.tb <- calc.0_rate(colon.uninfected.sim.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)
ileum.uninfected.sim.bin30.0rate.tb <- calc.0_rate(ileum.uninfected.sim.bin30.TL$matrices) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime)

hBECs.0_rate.per.bin <- rbind(cbind(data = "infected", hBECs.bin30.0rate.tb),
                              cbind(data = "simulated", hBECs.uninfected.sim.bin30.0rate.tb)) %>%
  ggplot(aes(y = value, x = pseudotime %>% factor(levels = 1:30), fill = data)) + geom_boxplot(outlier.shape = NA) + ylim(0, 1) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) + ggtitle("hBECs")

colon.0_rate.per.bin <- rbind(cbind(data = "infected", colon.bin30.0rate.tb),
                              cbind(data = "simulated", colon.uninfected.sim.bin30.0rate.tb)) %>%
  ggplot(aes(y = value, x = pseudotime %>% factor(levels = 1:30), fill = data)) + geom_boxplot(outlier.shape = NA) + ylim(0, 1) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) + ggtitle("colon")

ileum.0_rate.per.bin <- rbind(cbind(data = "infected", ileum.bin30.0rate.tb),
                              cbind(data = "simulated", ileum.uninfected.sim.bin30.0rate.tb)) %>%
  ggplot(aes(y = value, x = pseudotime %>% factor(levels = 1:30), fill = data)) + geom_boxplot(outlier.shape = NA) + ylim(0, 1) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) + ggtitle("ileum")

gg.0_rate.per.bin <- ggarrange(hBECs.0_rate.per.bin, colon.0_rate.per.bin, ileum.0_rate.per.bin,
                               nrow = 3, common.legend = TRUE, legend = "right") %>%
  annotate_figure(left = "proportion of zeros", bottom = "pseudotime bin", fig.lab = "a")

# clustering analysis with Seurat to ensure simulated cells cluster with uninfected cells

# hBECs

hBECs.uninfected.seurat <- CreateSeuratObject(hBECs.mtx$matrices$uninfected_0dpi[colnames(hBECs.mtx$matrices$uninfected_0dpi) %in%
                                                                                   colnames(hBECs.mtx.uninfected.down.sample)], project = "uninfected")
hBECs.sim.seurat <- CreateSeuratObject(hBECs.mtx.uninfected.down.sample, project = "sim")
hBECs.seurat <- merge(hBECs.uninfected.seurat, y = hBECs.sim.seurat, add.cell.ids = c("uninfected", "sim"))

hBECs.seurat$data <- c(rep("uninfected", ncol(hBECs.mtx.uninfected.down.sample)),
                       rep("simulated", ncol(hBECs.mtx.uninfected.down.sample)))
Idents(hBECs.seurat) <- "data"

hBECs.seurat <- NormalizeData(hBECs.seurat)
hBECs.seurat <- FindVariableFeatures(hBECs.seurat)
hBECs.seurat <- ScaleData(hBECs.seurat)
hBECs.seurat <- RunPCA(hBECs.seurat)
hBECs.seurat <- FindNeighbors(hBECs.seurat, dims = 1:10)
hBECs.seurat <- RunUMAP(hBECs.seurat, dims = 1:10)
hBECs.gg.umap <- DimPlot(hBECs.seurat, shuffle = TRUE)

# colon

colon.uninfected.seurat <- CreateSeuratObject(colon.mtx$matrices$uninfected_0hpi, project = "uninfected")
colon.sim.seurat <- CreateSeuratObject(colon.mtx.uninfected.down.sample, project = "sim")
colon.seurat <- merge(colon.uninfected.seurat, y = colon.sim.seurat, add.cell.ids = c("uninfected", "sim"))

colon.seurat$data <- c(rep("uninfected", ncol(colon.mtx$matrices$uninfected_0hpi)),
                       rep("simulated", ncol(colon.mtx.uninfected.down.sample)))
Idents(colon.seurat) <- "data"

colon.seurat <- NormalizeData(colon.seurat)
colon.seurat <- FindVariableFeatures(colon.seurat)
colon.seurat <- ScaleData(colon.seurat)
colon.seurat <- RunPCA(colon.seurat)
colon.seurat <- FindNeighbors(colon.seurat, dims = 1:10)
colon.seurat <- RunUMAP(colon.seurat, dims = 1:10)
colon.gg.umap <- DimPlot(colon.seurat, shuffle = TRUE)

# ileum

ileum.uninfected.seurat <- CreateSeuratObject(ileum.mtx$matrices$uninfected_0hpi, project = "uninfected")
ileum.sim.seurat <- CreateSeuratObject(ileum.mtx.uninfected.down.sample, project = "sim")
ileum.seurat <- merge(ileum.uninfected.seurat, y = ileum.sim.seurat, add.cell.ids = c("uninfected", "sim"))

ileum.seurat$data <- c(rep("uninfected", ncol(ileum.mtx$matrices$uninfected_0hpi)),
                       rep("simulated", ncol(ileum.mtx.uninfected.down.sample)))
Idents(ileum.seurat) <- "data"

ileum.seurat <- NormalizeData(ileum.seurat)
ileum.seurat <- FindVariableFeatures(ileum.seurat)
ileum.seurat <- ScaleData(ileum.seurat)
ileum.seurat <- RunPCA(ileum.seurat)
ileum.seurat <- FindNeighbors(ileum.seurat, dims = 1:10)
ileum.seurat <- RunUMAP(ileum.seurat, dims = 1:10)
ileum.gg.umap <- DimPlot(ileum.seurat, shuffle = TRUE)

# plot

gg.umap <- ggarrange(hBECs.gg.umap + ggtitle("hBECs") + theme(axis.title = element_blank()),
                     colon.gg.umap + ggtitle("colon") + theme(axis.title = element_blank()),
                     ileum.gg.umap + ggtitle("ileum") + theme(axis.title = element_blank()),
                     common.legend = TRUE, legend = "right", nrow = 1) %>%
  annotate_figure(left = "UMAP_2", bottom = "UMAP_1", fig.lab = "b")

gg.simulated <- ggarrange(gg.0_rate.per.bin, gg.umap, ncol = 1, heights = c(2, 1))

ggsave("figs_out/final/simulated_dataset_properties.pdf", gg.simulated, width = 6.85, height = 7)

# plot power laws for each pseudotime bin

hBECs.uninfected.sim.llmsd.bin30 <- setNames(hBECs.uninfected.sim.bin30.TL$log_log.mean_sd, nm = 1:length(hBECs.uninfected.sim.bin30.TL$matrices)) %>% bind_rows(.id = "bin")
colon.uninfected.sim.llmsd.bin30 <- setNames(colon.uninfected.sim.bin30.TL$log_log.mean_sd, nm = 1:length(colon.uninfected.sim.bin30.TL$matrices)) %>% bind_rows(.id = "bin")
ileum.uninfected.sim.llmsd.bin30 <- setNames(ileum.uninfected.sim.bin30.TL$log_log.mean_sd, nm = 1:length(ileum.uninfected.sim.bin30.TL$matrices)) %>% bind_rows(.id = "bin")

hBECs.uninfected.sim.llmsd.bin30$bin <- as.numeric(hBECs.uninfected.sim.llmsd.bin30$bin)
colon.uninfected.sim.llmsd.bin30$bin <- as.numeric(colon.uninfected.sim.llmsd.bin30$bin)
ileum.uninfected.sim.llmsd.bin30$bin <- as.numeric(ileum.uninfected.sim.llmsd.bin30$bin)

hBECs.uninfected.sim.llmsd.bin30$gene <- rownames(hBECs.uninfected.sim.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)
colon.uninfected.sim.llmsd.bin30$gene <- rownames(colon.uninfected.sim.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)
ileum.uninfected.sim.llmsd.bin30$gene <- rownames(ileum.uninfected.sim.llmsd.bin30) %>% gsub("\\.\\.\\..*", "", .)

hBECs.uninfected.sim.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", hBECs.uninfected.sim.llmsd.bin30$gene), "MT", "nuclear")
colon.uninfected.sim.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", colon.uninfected.sim.llmsd.bin30$gene), "MT", "nuclear")
ileum.uninfected.sim.llmsd.bin30$is.mt <- ifelse(grepl("^MT-", ileum.uninfected.sim.llmsd.bin30$gene), "MT", "nuclear")

hBECs.uninfected.sim.llmsd.bin30.gg <- ggplot(hBECs.uninfected.sim.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))
colon.uninfected.sim.llmsd.bin30.gg <- ggplot(colon.uninfected.sim.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("colon") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))
ileum.uninfected.sim.llmsd.bin30.gg <- ggplot(ileum.uninfected.sim.llmsd.bin30, aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("ileum") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))

# ANOVA

V.uninfected.sim.bin30.aov <- aov(lm(V ~ pseudotime_bin * cell + sparsity + n.genes, data = pseudotime_uninfected.sim.bin30.TL_params))
beta.uninfected.sim.bin30.aov <- aov(lm(beta ~ pseudotime_bin * cell + sparsity + n.genes, data = pseudotime_uninfected.sim.bin30.TL_params))

V.uninfected.sim.bin30.pes <- V.uninfected.sim.bin30.aov %>% effectsize::eta_squared()
beta.uninfected.sim.bin30.pes <- beta.uninfected.sim.bin30.aov %>% effectsize::eta_squared()

# segmented fit for simulated uninfected data

hBECs.uninfected.sim.bin30.threshold.TL <- pseudotime_bin(hBECs.mtx.uninfected.down.sample, zero.rate.threshold = .7, hBECs.n.bins)
colon.uninfected.sim.bin30.threshold.TL <- pseudotime_bin(colon.mtx.uninfected.down.sample, zero.rate.threshold = .7, colon.n.bins)
ileum.uninfected.sim.bin30.threshold.TL <- pseudotime_bin(ileum.mtx.uninfected.down.sample, zero.rate.threshold = .7, ileum.n.bins)

hBECs.uninfected.sim.bin30.threshold.seg.TL <- fit.seg.TL(hBECs.uninfected.sim.bin30.threshold.TL)
colon.uninfected.sim.bin30.threshold.seg.TL <- fit.seg.TL(colon.uninfected.sim.bin30.threshold.TL)
ileum.uninfected.sim.bin30.threshold.seg.TL <- fit.seg.TL(ileum.uninfected.sim.bin30.threshold.TL)

breakpoint.sparsity.uninfected.sim.bin30.tb <- rbind(hBECs.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "hBECs"),
                                                     colon.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "colon"),
                                                     ileum.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "ileum"))

colnames(breakpoint.sparsity.uninfected.sim.bin30.tb)[1] <- "pseudotime_bin"

# check plots

hBECs.llmsd.uninfected.sim.bin30 <- setNames(hBECs.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd, nm = 1:length(hBECs.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd)) %>% bind_rows(.id = "bin")
colon.llmsd.uninfected.sim.bin30 <- setNames(colon.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd, nm = 1:length(colon.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd)) %>% bind_rows(.id = "bin")
ileum.llmsd.uninfected.sim.bin30 <- setNames(ileum.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd, nm = 1:length(ileum.uninfected.sim.bin30.threshold.seg.TL$log_log.mean_sd)) %>% bind_rows(.id = "bin")

hBECs.llmsd.uninfected.sim.bin30$bin <- as.numeric(hBECs.llmsd.uninfected.sim.bin30$bin)
colon.llmsd.uninfected.sim.bin30$bin <- as.numeric(colon.llmsd.uninfected.sim.bin30$bin)
ileum.llmsd.uninfected.sim.bin30$bin <- as.numeric(ileum.llmsd.uninfected.sim.bin30$bin)

hBECs.llmsd.uninfected.sim.bin30$gene <- rownames(hBECs.llmsd.uninfected.sim.bin30) %>% gsub("\\.*", "", .)
colon.llmsd.uninfected.sim.bin30$gene <- rownames(colon.llmsd.uninfected.sim.bin30) %>% gsub("\\.*", "", .)
ileum.llmsd.uninfected.sim.bin30$gene <- rownames(ileum.llmsd.uninfected.sim.bin30) %>% gsub("\\.*", "", .)

hBECs.llmsd.uninfected.sim.bin30.gg <- ggplot(hBECs.llmsd.uninfected.sim.bin30, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")
colon.llmsd.uninfected.sim.bin30.gg <- ggplot(colon.llmsd.uninfected.sim.bin30, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("colon") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")
ileum.llmsd.uninfected.sim.bin30.gg <- ggplot(ileum.llmsd.uninfected.sim.bin30, aes(x = mean, y = sd)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("ileum") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -4, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none")

# ANOVA and other tests

breakpoint.sparsity.uninfected.sim.bin30.tb <- rbind(hBECs.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "hBECs"),
                                                     colon.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "colon"),
                                                     ileum.uninfected.sim.bin30.threshold.seg.TL$params %>% cbind(cell = "ileum")) %>%
  separate(V, into = c("V1", "V2"), sep = ",") %>%
  separate(beta, into = c("beta1", "beta2"), sep = ",")

breakpoint.sparsity.uninfected.sim.bin30.tb$cell <- factor(breakpoint.sparsity.uninfected.sim.bin30.tb$cell, levels = c("hBECs", "colon", "ileum"))

colnames(breakpoint.sparsity.uninfected.sim.bin30.tb)[1] <- "pseudotime_bin"

breakpoint.sparsity.uninfected.sim.bin30.tb$pseudotime_bin <- breakpoint.sparsity.uninfected.sim.bin30.tb$pseudotime_bin %>% as.numeric()
breakpoint.sparsity.uninfected.sim.bin30.tb$breakpoint <- breakpoint.sparsity.uninfected.sim.bin30.tb$breakpoint %>% as.numeric()
breakpoint.sparsity.uninfected.sim.bin30.tb$V1 <- breakpoint.sparsity.uninfected.sim.bin30.tb$V1 %>% as.numeric()
breakpoint.sparsity.uninfected.sim.bin30.tb$V2 <- breakpoint.sparsity.uninfected.sim.bin30.tb$V2 %>% as.numeric()
breakpoint.sparsity.uninfected.sim.bin30.tb$beta1 <- breakpoint.sparsity.uninfected.sim.bin30.tb$beta1 %>% as.numeric()
breakpoint.sparsity.uninfected.sim.bin30.tb$beta2 <- breakpoint.sparsity.uninfected.sim.bin30.tb$beta2 %>% as.numeric()

# test breakpoint

uninfected.sim.breakpoint.lm <- lm(breakpoint ~ pseudotime_bin * cell + n.genes + sparsity,
                                   data = breakpoint.sparsity.uninfected.sim.bin30.tb)

uninfected.sim.breakpoint.aov <- aov(uninfected.sim.breakpoint.lm)

uninfected.sim.breakpoint.pes <- effectsize::eta_squared(uninfected.sim.breakpoint.lm)

# test estimated V's and beta's with the breakpoint as a covariate

uninfected.sim.V1.lm <- lm(V1 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
                           data = breakpoint.sparsity.uninfected.sim.bin30.tb)
uninfected.sim.beta1.lm <- lm(beta1 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
                              data = breakpoint.sparsity.uninfected.sim.bin30.tb)

uninfected.sim.V2.lm <- lm(V2 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
                           data = breakpoint.sparsity.uninfected.sim.bin30.tb)
uninfected.sim.beta2.lm <- lm(beta2 ~ pseudotime_bin * cell + breakpoint + n.genes + sparsity,
                              data = breakpoint.sparsity.uninfected.sim.bin30.tb)

uninfected.sim.V1.aov <- aov(uninfected.sim.V1.lm)
uninfected.sim.beta1.aov <- aov(uninfected.sim.beta1.lm)
uninfected.sim.V2.aov <- aov(uninfected.sim.V2.lm)
uninfected.sim.beta2.aov <- aov(uninfected.sim.beta2.lm)

uninfected.sim.V1.pes <- effectsize::eta_squared(uninfected.sim.V1.lm)
uninfected.sim.beta1.pes <- effectsize::eta_squared(uninfected.sim.beta1.lm)
uninfected.sim.V2.pes <- effectsize::eta_squared(uninfected.sim.V2.lm)
uninfected.sim.beta2.pes <- effectsize::eta_squared(uninfected.sim.beta2.lm)

# plot

breakpoint.sparsity.uninfected.sim.bin30.gg.tb <- breakpoint.sparsity.uninfected.sim.bin30.tb %>%
  dplyr::select(pseudotime_bin, V1, V2, beta1, beta2, breakpoint, cell) %>%
  pivot_longer(-c(pseudotime_bin, cell))

breakpoint.sparsity.uninfected.sim.bin30.gg.tb$name <- factor(breakpoint.sparsity.uninfected.sim.bin30.gg.tb$name, levels = c("breakpoint", "V1", "beta1", "V2", "beta2"))

uninfected.sim.npsi1.breakpoint.gg <- ggplot(breakpoint.sparsity.uninfected.sim.bin30.gg.tb, aes(pseudotime_bin, value)) + geom_line() +
  facet_grid(cols = vars(cell), rows = vars(name), scales = "free") + theme(axis.title.y = element_blank()) + xlab("pseudotime bin")

# save data

save(hBECs.mtx.uninfected.down.sample,
     colon.mtx.uninfected.down.sample,
     ileum.mtx.uninfected.down.sample,
     hBECs.uninfected.sim.bin30.TL,
     colon.uninfected.sim.bin30.TL,
     ileum.uninfected.sim.bin30.TL,
     file = "RData/TPT_MT/sim_cell_infection_pseudotime_bin_matrices_TL.RData")

# figure for article

# plot bins

breakpoint_bin30.npsi1 <- breakpoint.sparsity.bin30.tb %>% mutate(V = V1, beta = beta1, model = "segmented") %>%
  dplyr::select(pseudotime_bin, V, beta, R.squared, n.cells, n.genes, sparsity, model, cell)

pseudotime_bin30.TL_params.unseg.seg <- rbind(pseudotime_bin30.TL_params %>% mutate(model := "unsegmented"),
                                              breakpoint_bin30.npsi1)

pseudotime_bin30.TL_params.unseg.seg$shape <- ifelse(paste0(pseudotime_bin30.TL_params.unseg.seg$pseudotime_bin,
                                                            pseudotime_bin30.TL_params.unseg.seg$cell) %in%
                                                       c("1hBECs", "18hBECs", "23hBECs", "30hBECs",
                                                          "1colon", "13colon", "27colon", "30colon",
                                                          "1ileum", "21ileum", "28ileum", "30ileum"), "a", "b")

pseudotime_bin30.TL_params$shape <- ifelse(paste0(pseudotime_bin30.TL_params$pseudotime_bin,
                                                  pseudotime_bin30.TL_params$cell) %in%
                                             c("1hBECs", "18hBECs", "23hBECs", "30hBECs",
                                               "1colon", "13colon", "27colon", "30colon",
                                               "1ileum", "21ileum", "28ileum", "30ileum"), "a", "b")

# plot pseudotime_bin30.TL_params or pseudotime_bin30.TL_params.unseg.seg

bins.gg <- ggplot(pseudotime_bin30.TL_params %>%
                    dplyr::select(c(pseudotime_bin, V, beta, cell, shape)) %>%
                    pivot_longer(-c(cell, pseudotime_bin, shape)),
                  aes(x = pseudotime_bin, y = value, color = cell)) +
  geom_point(aes(shape = shape), size = 1.2) + scale_shape_manual(values = c(18, NA)) + guides(shape = "none") +
  facet_wrap(~name %>% factor(levels = c("V", "beta")), scales = "free", strip.position = "left", ncol = 1) +
  geom_path(linewidth = .25) +
  theme(axis.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12))

bins.sim.gg <- ggplot(pseudotime_uninfected.sim.bin30.TL_params %>%
                        dplyr::select(c(pseudotime_bin, V, beta, cell)) %>%
                        pivot_longer(-c(cell, pseudotime_bin)),
                      aes(x = pseudotime_bin, y = value, color = cell)) +
  facet_wrap(~name %>% factor(levels = c("V", "beta")), scales = "free", strip.position = "left", ncol = 1,
             labeller = c(V = "",
                          beta = "") %>% as_labeller()) +
  geom_path(linewidth = .25) +
  theme(axis.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12))

all.bins.gg <- rbind(cbind(pseudotime_bin30.TL_params, data = "infected"),
      cbind(pseudotime_uninfected.sim.bin30.TL_params, shape = "b", data = "simulated")) %>% 
  dplyr::select(c(pseudotime_bin, V, beta, cell, shape, data)) %>%
  pivot_longer(-c(cell, pseudotime_bin, shape, data)) %>%
  ggplot(aes(x = pseudotime_bin, y = value, color = cell)) +
  geom_path(linewidth = .25) +
  geom_point(aes(shape = shape), size = 1.2) +
  scale_shape_manual(values = c(19, NA)) +
  guides(shape = "none") +
  facet_grid(vars(name %>% factor(levels = c("V", "beta"))), vars(data), scales = "free", switch = "y",
             labeller = label_parsed) +
  xlab("pseudotime bin") +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  labs(tag = "a")

# power laws
# transform to log-log base 10!

hBECs.llmsd.bin30.gg.a <- ggplot(hBECs.llmsd.bin30 %>%
                                   mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)) %>% 
                                   dplyr::filter(bin %in% c(1, 18, 23, 30)),
                                 aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("hBECs") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -1.7, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("orange", "green", "black"))

colon.llmsd.bin30.gg.a <- ggplot(colon.llmsd.bin30 %>%
                                   mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)) %>%
                                   dplyr::filter(bin %in% c(1, 13, 27, 30)),
                                 aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("colon") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -1.7, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))

ileum.llmsd.bin30.gg.a <- ggplot(ileum.llmsd.bin30 %>%
                                   mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)) %>%
                                   dplyr::filter(bin %in% c(1, 21, 28, 30)),
                                 aes(x = mean, y = sd, color = is.mt)) + xlab("mean (log)") + ylab("standard deviation (log)") + ggtitle("ileum") +
  geom_point(size = .01) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = .25) + geom_abline(intercept = -1.7, slope = .5, color = "blue", linewidth = .25) +
  facet_wrap(~bin, ncol = 5) + theme(legend.position = "none") + scale_color_manual(values = c("green", "black"))

selected.bins.gg <- ggarrange(hBECs.llmsd.bin30.gg.a + theme(axis.title = element_blank(),
                                                             title = element_text(size = 7),
                                                             axis.text = element_text(size = 5)),
                              colon.llmsd.bin30.gg.a + theme(axis.title = element_blank(),
                                                             title = element_text(size = 7),
                                                             axis.text = element_text(size = 5)),
                              ileum.llmsd.bin30.gg.a + theme(axis.title = element_blank(),
                                                             title = element_text(size = 7),
                                                             axis.text = element_text(size = 5)), ncol = 1) %>%
  annotate_figure(left = "standard deviation (log)", bottom = "mean (log)", fig.lab = " b", fig.lab.size = 14)

gg2 <- ggarrange(all.bins.gg,
                 selected.bins.gg,
                 ncol = 1, heights = c(.4, .6))

ggsave("figs_out/final/selected_bins.pdf", gg2, width = 6.85, height = 9.21)
