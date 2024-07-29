source("source_me.R")

# read/create rank matrices

load("RData/TPT_MT/cell_infection_pseudotime_bin_matrices_TL.RData")
load("RData/TPT_MT/cell_infection_pseudotime_bin_matrices_seg_TL.RData")

load(paste0("RData/", out.dir, "cell_infection_pseudotime_matrices_TL.RData"))
load("RData/TPT_MT/sim_cell_infection_pseudotime_bin_matrices_TL.RData")

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

# rank 

hBECs.rank <- process.rank(matrices = hBECs.mtx$matrices[c("infected", "uninfected_0dpi")], method = "random")
colon.rank <- process.rank(matrices = colon.mtx$matrices[c("infected", "uninfected_0hpi")], method = "random")
ileum.rank <- process.rank(matrices = ileum.mtx$matrices[c("infected", "uninfected_0hpi")], method = "random")

# RSI of the mean gene expression at each pseudotime point

# create matrices with the gene mean for each pseudotime bin

hBECs.bin30.no_filter.TL <- pseudotime_bin(hBECs.mtx$matrices$infected, zero.rate.threshold = NULL, n.bins = hBECs.n.bins)
colon.bin30.no_filter.TL <- pseudotime_bin(colon.mtx$matrices$infected, zero.rate.threshold = NULL, n.bins = colon.n.bins)
ileum.bin30.no_filter.TL <- pseudotime_bin(ileum.mtx$matrices$infected, zero.rate.threshold = NULL, n.bins = ileum.n.bins)

hBECs.bin30_mean_sd.mtx <- hBECs.bin30.no_filter.TL %>% get.mean.sd()
colon.bin30_mean_sd.mtx <- colon.bin30.no_filter.TL %>% get.mean.sd()
ileum.bin30_mean_sd.mtx <- ileum.bin30.no_filter.TL %>% get.mean.sd()

hBECs.bin30_mean.rank <- process.rank(matrices = list(infected = hBECs.bin30_mean_sd.mtx$mean), method = "random")
colon.bin30_mean.rank <- process.rank(matrices = list(infected = colon.bin30_mean_sd.mtx$mean), method = "random")
ileum.bin30_mean.rank <- process.rank(matrices = list(infected = ileum.bin30_mean_sd.mtx$mean), method = "random")

hBECs.bin30_mean.RSI <- calc.RSI(hBECs.bin30_mean.rank)
colon.bin30_mean.RSI <- calc.RSI(colon.bin30_mean.rank)
ileum.bin30_mean.RSI <- calc.RSI(ileum.bin30_mean.rank)

# RSI bootstrap, only analyse genes with 0 rate = 0

hBECs.bin30_mean.0rate <- calc.0_rate(list(infected = hBECs.bin30_mean_sd.mtx$mean))
colon.bin30_mean.0rate <- calc.0_rate(list(infected = colon.bin30_mean_sd.mtx$mean))
ileum.bin30_mean.0rate <- calc.0_rate(list(infected = ileum.bin30_mean_sd.mtx$mean))

hBECs.bin30_mean.RSI.boot <- calc.RSI(hBECs.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(hBECs.bin30_mean.0rate$infected[hBECs.bin30_mean.0rate$infected == 0]))
colon.bin30_mean.RSI.boot <- calc.RSI(colon.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(colon.bin30_mean.0rate$infected[colon.bin30_mean.0rate$infected == 0]))
ileum.bin30_mean.RSI.boot <- calc.RSI(ileum.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(ileum.bin30_mean.0rate$infected[ileum.bin30_mean.0rate$infected == 0]))

hBECs.bin30_mean.RSI.boot.temporal_genes <- hBECs.bin30_mean.RSI.boot$infected[hBECs.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,]
colon.bin30_mean.RSI.boot.temporal_genes <- colon.bin30_mean.RSI.boot$infected[colon.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,]
ileum.bin30_mean.RSI.boot.temporal_genes <- ileum.bin30_mean.RSI.boot$infected[ileum.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,]

# Taylor's plot with the means of each bin

means.TL <- fit.TL(list(hBECs = hBECs.bin30_mean_sd.mtx$mean,
                        colon = colon.bin30_mean_sd.mtx$mean,
                        ileum = ileum.bin30_mean_sd.mtx$mean))

plot.power_law(means.TL$log_log.mean_sd, list(hBECs = setNames(ifelse(hBECs.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,
                                                                      "< 0.05", ">= 0.05"),
                                                               nm = hBECs.bin30_mean.RSI.boot$infected %>% rownames()),
                                              colon = setNames(ifelse(colon.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,
                                                                      "< 0.05", ">= 0.05"),
                                                               nm = colon.bin30_mean.RSI.boot$infected %>% rownames()),
                                              ileum = setNames(ifelse(ileum.bin30_mean.RSI.boot$infected$p.adjust.survival < .05,
                                                                      "< 0.05", ">= 0.05"),
                                                               nm = ileum.bin30_mean.RSI.boot$infected %>% rownames())),
               "FDR", metadata.range = NULL) + theme(legend.position = "none") +
  scale_color_manual(values = c("tomato", "grey60"))

# random matrices

set.seed(1024)
hBECs.random.mtx <- matrix(runif(ncol(hBECs.rank$infected) * nrow(hBECs.rank$infected), min = -1, max = 1),
                           nrow = nrow(hBECs.rank$infected)) %>% as.data.frame()
rownames(hBECs.random.mtx) <- rownames(hBECs.rank$infected)

set.seed(1024)
colon.random.mtx <- matrix(runif(ncol(colon.rank$infected) * nrow(colon.rank$infected), min = -1, max = 1),
                           nrow = nrow(colon.rank$infected)) %>% as.data.frame()
rownames(colon.random.mtx) <- rownames(colon.rank$infected)

set.seed(1024)
ileum.random.mtx <- matrix(runif(ncol(ileum.rank$infected) * nrow(ileum.rank$infected), min = -1, max = 1),
                           nrow = nrow(ileum.rank$infected)) %>% as.data.frame()
rownames(ileum.random.mtx) <- rownames(ileum.rank$infected)

# RSI boot random test

#hBECs.RSI.rand <- process.rank(matrices = list(rand = hBECs.random.mtx), method = "random") %>% calc.RSI(shuffle.rep = 1000)

# also test RSI for uninfected cells as a null hypothesis

hBECs.bin30.null.TL <- pseudotime_bin(hBECs.mtx$matrices$uninfected_0dpi, zero.rate.threshold = NULL)
colon.bin30.null.TL <- pseudotime_bin(colon.mtx$matrices$uninfected_0hpi, zero.rate.threshold = NULL)
ileum.bin30.null.TL <- pseudotime_bin(ileum.mtx$matrices$uninfected_0hpi, zero.rate.threshold = NULL)

hBECs.bin30.null_mean_sd.mtx <- hBECs.bin30.null.TL %>% get.mean.sd()
colon.bin30.null_mean_sd.mtx <- colon.bin30.null.TL %>% get.mean.sd()
ileum.bin30.null_mean_sd.mtx <- ileum.bin30.null.TL %>% get.mean.sd()

hBECs.bin30.null_mean.rank <- process.rank(matrices = list(uninfected = hBECs.bin30.null_mean_sd.mtx$mean), method = "random")
colon.bin30.null_mean.rank <- process.rank(matrices = list(uninfected = colon.bin30.null_mean_sd.mtx$mean), method = "random")
ileum.bin30.null_mean.rank <- process.rank(matrices = list(uninfected = ileum.bin30.null_mean_sd.mtx$mean), method = "random")

hBECs.bin30.null_mean.0rate <- calc.0_rate(list(uninfected = hBECs.bin30.null_mean_sd.mtx$mean))
colon.bin30.null_mean.0rate <- calc.0_rate(list(uninfected = colon.bin30.null_mean_sd.mtx$mean))
ileum.bin30.null_mean.0rate <- calc.0_rate(list(uninfected = ileum.bin30.null_mean_sd.mtx$mean))

hBECs.bin30.null_mean.RSI.boot <- calc.RSI(hBECs.bin30.null_mean.rank, shuffle.rep = 1000, genes.to.keep = names(hBECs.bin30.null_mean.0rate$uninfected[hBECs.bin30.null_mean.0rate$uninfected == 0]))
colon.bin30.null_mean.RSI.boot <- calc.RSI(colon.bin30.null_mean.rank, shuffle.rep = 1000, genes.to.keep = names(colon.bin30.null_mean.0rate$uninfected[colon.bin30.null_mean.0rate$uninfected == 0]))
ileum.bin30.null_mean.RSI.boot <- calc.RSI(ileum.bin30.null_mean.rank, shuffle.rep = 1000, genes.to.keep = names(ileum.bin30.null_mean.0rate$uninfected[ileum.bin30.null_mean.0rate$uninfected == 0]))

hBECs.bin30.null_mean.RSI.boot.temporal_genes <- hBECs.bin30.null_mean.RSI.boot$uninfected[hBECs.bin30.null_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                            !is.nan(hBECs.bin30.null_mean.RSI.boot$uninfected$p.adjust),]
colon.bin30.null_mean.RSI.boot.temporal_genes <- colon.bin30.null_mean.RSI.boot$uninfected[colon.bin30.null_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                             !is.nan(colon.bin30.null_mean.RSI.boot$uninfected$p.adjust),]
ileum.bin30.null_mean.RSI.boot.temporal_genes <- ileum.bin30.null_mean.RSI.boot$uninfected[ileum.bin30.null_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                             !is.nan(ileum.bin30.null_mean.RSI.boot$uninfected$p.adjust),]

RSI.boot.bin30.null_intersect <- Reduce(intersect, list(hBECs.bin30.null_mean.RSI.boot.temporal_genes %>% rownames(),
                                                        colon.bin30.null_mean.RSI.boot.temporal_genes %>% rownames(),
                                                        ileum.bin30.null_mean.RSI.boot.temporal_genes %>% rownames())) %>% sort()

# simulated infected dataset

hBECs.uninfected.sim.bin30_mean_sd.mtx <- hBECs.uninfected.sim.bin30.TL %>% get.mean.sd()
colon.uninfected.sim.bin30_mean_sd.mtx <- colon.uninfected.sim.bin30.TL %>% get.mean.sd()
ileum.uninfected.sim.bin30_mean_sd.mtx <- ileum.uninfected.sim.bin30.TL %>% get.mean.sd()

hBECs.uninfected.sim.bin30_mean.rank <- process.rank(matrices = list(uninfected = hBECs.uninfected.sim.bin30_mean_sd.mtx$mean), method = "random")
colon.uninfected.sim.bin30_mean.rank <- process.rank(matrices = list(uninfected = colon.uninfected.sim.bin30_mean_sd.mtx$mean), method = "random")
ileum.uninfected.sim.bin30_mean.rank <- process.rank(matrices = list(uninfected = ileum.uninfected.sim.bin30_mean_sd.mtx$mean), method = "random")

hBECs.uninfected.sim.bin30_mean.0rate <- calc.0_rate(list(uninfected = hBECs.uninfected.sim.bin30_mean_sd.mtx$mean))
colon.uninfected.sim.bin30_mean.0rate <- calc.0_rate(list(uninfected = colon.uninfected.sim.bin30_mean_sd.mtx$mean))
ileum.uninfected.sim.bin30_mean.0rate <- calc.0_rate(list(uninfected = ileum.uninfected.sim.bin30_mean_sd.mtx$mean))

hBECs.uninfected.sim.bin30_mean.RSI.boot <- calc.RSI(hBECs.uninfected.sim.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(hBECs.uninfected.sim.bin30_mean.0rate$uninfected[hBECs.uninfected.sim.bin30_mean.0rate$uninfected == 0]))
colon.uninfected.sim.bin30_mean.RSI.boot <- calc.RSI(colon.uninfected.sim.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(colon.uninfected.sim.bin30_mean.0rate$uninfected[colon.uninfected.sim.bin30_mean.0rate$uninfected == 0]))
ileum.uninfected.sim.bin30_mean.RSI.boot <- calc.RSI(ileum.uninfected.sim.bin30_mean.rank, shuffle.rep = 1000, genes.to.keep = names(ileum.uninfected.sim.bin30_mean.0rate$uninfected[ileum.uninfected.sim.bin30_mean.0rate$uninfected == 0]))

hBECs.uninfected.sim.bin30_mean.RSI.boot.temporal_genes <- hBECs.uninfected.sim.bin30_mean.RSI.boot$uninfected[hBECs.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                                                 !is.nan(hBECs.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust),]
colon.uninfected.sim.bin30_mean.RSI.boot.temporal_genes <- colon.uninfected.sim.bin30_mean.RSI.boot$uninfected[colon.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                                                 !is.nan(colon.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust),]
ileum.uninfected.sim.bin30_mean.RSI.boot.temporal_genes <- ileum.uninfected.sim.bin30_mean.RSI.boot$uninfected[ileum.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust.survival < .05 &
                                                                                                                 !is.nan(ileum.uninfected.sim.bin30_mean.RSI.boot$uninfected$p.adjust),]

RSI.boot.uninfected.sim.bin30_intersect <- Reduce(intersect, list(hBECs.uninfected.sim.bin30_mean.RSI.boot.temporal_genes %>% rownames(),
                                                                  colon.uninfected.sim.bin30_mean.RSI.boot.temporal_genes %>% rownames(),
                                                                  ileum.uninfected.sim.bin30_mean.RSI.boot.temporal_genes %>% rownames())) %>% sort()

# similar genes

RSI.boot.bin30_intersect <- Reduce(intersect, list(hBECs.bin30_mean.RSI.boot.temporal_genes %>% rownames(),
                                                   colon.bin30_mean.RSI.boot.temporal_genes %>% rownames(),
                                                   ileum.bin30_mean.RSI.boot.temporal_genes %>% rownames())) %>% sort()

# get their RSI

RSI.boot.bin30_intersect.RSI.values <- lapply(list(hBECs = hBECs.bin30_mean.RSI.boot.temporal_genes,
                                                   colon = colon.bin30_mean.RSI.boot.temporal_genes,
                                                   ileum = ileum.bin30_mean.RSI.boot.temporal_genes),
                                              function(x) setNames(x[rownames(x) %in% RSI.boot.bin30_intersect, "RSI"],
                                                                   nm = x[rownames(x) %in% RSI.boot.bin30_intersect,] %>% rownames())[order(x[rownames(x) %in% RSI.boot.bin30_intersect, "RSI"])])

# plot RSI, PSI (RSIboot) and FDR

RSI.RSIboot.FDR.tb <- rbind(cbind(cell = "hBECs", hBECs.bin30_mean.RSI.boot$infected),
                            cbind(cell = "colon", colon.bin30_mean.RSI.boot$infected),
                            cbind(cell = "ileum", ileum.bin30_mean.RSI.boot$infected))

RSI.RSIboot.FDR.gg <- RSI.RSIboot.FDR.tb %>%
  ggplot(aes(x = RSI, y = RSI.boot)) +
  facet_wrap(~cell %>% factor(levels = c("hBECs", "colon", "ileum"))) +
  geom_point(aes(color = ifelse(p.adjust.survival < .05, "a", "b"),
                 shape = ifelse(p.adjust.survival < .05, "a", "b")),
             size = .75, shape = 1, stroke = .08) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("tomato", "grey60")) +
  
  ylab("PSI") + labs(tag = "a")

# plot GO

RSI.boot.bin30_intersect.GO <- enrichGO(RSI.boot.bin30_intersect,
                                        org.Hs.eg.db, keyType = "SYMBOL",
                                        ont = "BP") %>% as.data.frame()

RSI.boot.bin30_intersect.bar.GO <- gg.bar.GO(RSI.boot.bin30_intersect.GO[1:min(50, nrow(RSI.boot.bin30_intersect.GO)),],
                                             gene.values = RSI.boot.bin30_intersect.RSI.values, legend = "median\nRSI") %>%
  annotate_figure(fig.lab = "  b", fig.lab.size = 14)

# save plot

RSI.gg <- ggarrange(RSI.RSIboot.FDR.gg + theme(aspect.ratio = 1,
                                               plot.tag.position = c(-.109, 1)),
                    RSI.boot.bin30_intersect.bar.GO,
                    heights = c(1, 2.8),
                    ncol = 1)

ggsave("figs_out/final/RSI_PSI_GO.eps", device = "eps",  RSI.gg, width = 6.85, height = 9)

# save RSI and GO results

get.RSI.info <- function(x) x %>% dplyr::select(RSI, RSI.rep.mean, p.value.survival, p.adjust.survival) %>%
  mutate(PSI = RSI / RSI.rep.mean, gene = rownames(x))

write_xlsx(list(hBECs_RSI = hBECs.bin30_mean.RSI.boot$infected %>% get.RSI.info(),
                colon_RSI = colon.bin30_mean.RSI.boot$infected %>% get.RSI.info(),
                ileum_RSI = ileum.bin30_mean.RSI.boot$infected %>% get.RSI.info(),
                
                hBECs_simulated_RSI = hBECs.uninfected.sim.bin30_mean.RSI.boot$uninfected %>% get.RSI.info(),
                colon_simulated_RSI = colon.uninfected.sim.bin30_mean.RSI.boot$uninfected %>% get.RSI.info(),
                ileum_simulated_RSI = ileum.uninfected.sim.bin30_mean.RSI.boot$uninfected %>% get.RSI.info(),
                
                hBECs_uninfected_RSI = hBECs.bin30.null_mean.RSI.boot$uninfected %>% get.RSI.info(),
                colon_uninfected_RSI = colon.bin30.null_mean.RSI.boot$uninfected %>% get.RSI.info(),
                ileum_uninfected_RSI = ileum.bin30.null_mean.RSI.boot$uninfected %>% get.RSI.info(),
                
                GO_results = RSI.boot.bin30_intersect.GO %>%
                  add.value.to.GO(RSI.boot.bin30_intersect.RSI.values)), "RSI_supplementary.xlsx")

# calculate Hurst exponent

# create more control matrices

set.seed(1024)
hBECs.random.mtx <- matrix(runif(ncol(hBECs.rank$infected) * nrow(hBECs.rank$infected), min = -1, max = 1),
                           nrow = nrow(hBECs.rank$infected)) %>% as.data.frame()
rownames(hBECs.random.mtx) <- rownames(hBECs.rank$infected)
colnames(hBECs.random.mtx) <- colnames(hBECs.rank$infected)
hBECs.sim.random.rank <- process.rank(matrices = list(simulated = hBECs.mtx.uninfected.down.sample,
                                                      random = hBECs.random.mtx))
set.seed(1024)
hBECs.infected.shuffled <- list(shuffled = as.data.frame(hBECs.rank$infected)[sample(1:ncol(hBECs.rank$infected))])

set.seed(1024)
colon.random.mtx <- matrix(runif(ncol(colon.rank$infected) * nrow(colon.rank$infected), min = -1, max = 1),
                           nrow = nrow(colon.rank$infected)) %>% as.data.frame()
rownames(colon.random.mtx) <- rownames(colon.rank$infected)
colnames(colon.random.mtx) <- colnames(colon.rank$infected)
colon.sim.random.rank <- process.rank(matrices = list(simulated = colon.mtx.uninfected.down.sample,
                                                      random = colon.random.mtx))
set.seed(1024)
colon.infected.shuffled <- list(shuffled = as.data.frame(colon.rank$infected)[sample(1:ncol(colon.rank$infected))])

set.seed(1024)
ileum.random.mtx <- matrix(runif(ncol(ileum.rank$infected) * nrow(ileum.rank$infected), min = -1, max = 1),
                           nrow = nrow(ileum.rank$infected)) %>% as.data.frame()
rownames(ileum.random.mtx) <- rownames(ileum.rank$infected)
colnames(ileum.random.mtx) <- colnames(ileum.rank$infected)
ileum.sim.random.rank <- process.rank(matrices = list(simulated = ileum.mtx.uninfected.down.sample,
                                                      random = ileum.random.mtx))
set.seed(1024)
ileum.infected.shuffled <- list(shuffled = as.data.frame(ileum.rank$infected)[sample(1:ncol(ileum.rank$infected))])

# calculate exponent for RSI boot genes

set.seed(1024)
hBECs.hurst.rank.RSI.boot <- fit.hurst(c(hBECs.rank,
                                         hBECs.infected.shuffled,
                                         hBECs.sim.random.rank),
                                       hBECs.bin30_mean.RSI.boot.temporal_genes %>% rownames())

set.seed(1024)
colon.hurst.rank.RSI.boot <- fit.hurst(c(colon.rank,
                                         colon.infected.shuffled,
                                         colon.sim.random.rank),
                                       colon.bin30_mean.RSI.boot.temporal_genes %>% rownames())

set.seed(1024)
ileum.hurst.rank.RSI.boot <- fit.hurst(c(ileum.rank,
                                         ileum.infected.shuffled,
                                         ileum.sim.random.rank),
                                       ileum.bin30_mean.RSI.boot.temporal_genes %>% rownames())

# plot

hBECs.hurst.rank.gg <- hBECs.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
  ggplot(aes(x = value, color = data)) + geom_density(alpha = .2, linewidth = 1) + facet_wrap(~name) + 
  ggtitle("hBECs") + theme(axis.title.x = element_blank())

colon.hurst.rank.gg <- colon.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
  ggplot(aes(x = value, color = data)) + geom_density(alpha = .2) + facet_wrap(~name) +
  ggtitle("colon") + theme(axis.title.x = element_blank())

ileum.hurst.rank.gg <- ileum.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
  ggplot(aes(x = value, color = data)) + geom_density(alpha = .2) + facet_wrap(~name) +
  ggtitle("ileum") + theme(axis.title.x = element_blank())

hurst.all_data.gg <- rbind(hBECs.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
                             cbind(cell = "hBECs"),
                           colon.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
                             cbind(cell = "colon"),
                           ileum.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% dplyr::select(-c(Ht, Hrs, Hs, gene)) %>% pivot_longer(-data) %>%
                             cbind(cell = "ileum")) %>% mutate(data := data %>% gsub("_0.*", "", .)) %>%
  ggplot(aes(x = value, color = data)) + geom_density(linewidth = .7) + facet_wrap(~factor(cell, levels = c("hBECs", "colon", "ileum"))) + 
  xlab("H")# + labs(tag = "a")

# calculate H for significant RSI boot genes vs all genes for each cell type

hBECs.hurst.rank.non_RSI.boot <- fit.hurst(c(hBECs.rank["infected"],
                                             hBECs.sim.random.rank["simulated"]),
                                           rownames(hBECs.rank$infected)[!rownames(hBECs.rank$infected) %in%
                                                                           rownames(hBECs.bin30_mean.RSI.boot.temporal_genes)])
colon.hurst.rank.non_RSI.boot <- fit.hurst(c(colon.rank["infected"],
                                             colon.sim.random.rank["simulated"]),
                                           rownames(colon.rank$infected)[!rownames(colon.rank$infected) %in%
                                                                           rownames(colon.bin30_mean.RSI.boot.temporal_genes)])
ileum.hurst.rank.non_RSI.boot <- fit.hurst(c(ileum.rank["infected"],
                                             ileum.sim.random.rank["simulated"]),
                                           rownames(ileum.rank$infected)[!rownames(ileum.rank$infected) %in%
                                                                           rownames(ileum.bin30_mean.RSI.boot.temporal_genes)])

hBECs.hurst.rank <- rbind(hBECs.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = "<0.05"),
                          hBECs.hurst.rank.non_RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = ">=0.05"))

colon.hurst.rank <- rbind(colon.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = "<0.05"),
                          colon.hurst.rank.non_RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = ">=0.05"))

ileum.hurst.rank <- rbind(ileum.hurst.rank.RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = "<0.05"),
                          ileum.hurst.rank.non_RSI.boot %>% bind_rows(.id = "gene") %>% cbind(`RSIboot FDR` = ">=0.05"))

hBECs.hurst.rank.RSI_vs_non.gg <- hBECs.hurst.rank %>% filter(data == "infected") %>% 
  dplyr::select(-c(Ht, Hrs, Hs, gene, data)) %>% pivot_longer(-`RSIboot FDR`) %>%
  ggplot(aes(x = value, fill = `RSIboot FDR`, color = `RSIboot FDR`)) + geom_density(alpha = .3) + facet_wrap(~name) +
  ggtitle("hBECs") + theme(axis.title.x = element_blank())

colon.hurst.rank.RSI_vs_non.gg <- colon.hurst.rank %>% filter(data == "infected") %>% 
  dplyr::select(-c(Ht, Hrs, Hs, gene, data)) %>% pivot_longer(-`RSIboot FDR`) %>%
  ggplot(aes(x = value, fill = `RSIboot FDR`, color = `RSIboot FDR`)) + geom_density(alpha = .3) + facet_wrap(~name) +
  ggtitle("colon") + theme(axis.title.x = element_blank())

ileum.hurst.rank.RSI_vs_non.gg <- ileum.hurst.rank %>% filter(data == "infected") %>% 
  dplyr::select(-c(Ht, Hrs, Hs, gene, data)) %>% pivot_longer(-`RSIboot FDR`) %>%
  ggplot(aes(x = value, fill = `RSIboot FDR`, color = `RSIboot FDR`)) + geom_density(alpha = .3) + facet_wrap(~name) +
  ggtitle("ileum") + theme(axis.title.x = element_blank())

# plot all genes

hurst.rank.all.gg <- rbind(hBECs.hurst.rank %>% filter(data == "infected") %>% cbind(cell = "hBECs"),
                           colon.hurst.rank %>% filter(data == "infected") %>% cbind(cell = "colon"),
                           ileum.hurst.rank %>% filter(data == "infected") %>% cbind(cell = "ileum")) %>% 
  dplyr::select(-c(Ht, Hrs, Hs, gene, data)) %>% pivot_longer(-c(cell, `RSIboot FDR`)) %>%
  ggplot(aes(x = value, fill = cell, color = cell)) + geom_density(alpha = .3) + facet_wrap(~name) +
  theme(axis.title.x = element_blank()) +
  scale_fill_discrete(breaks = c("hBECs", "colon", "ileum")) + scale_color_discrete(breaks = c("hBECs", "colon", "ileum"))

# pairwise gene comparison, infected vs simulated

hurst.infected.vs.simulated.tb <- rbind(hBECs.hurst.rank %>% cbind(cell = "hBECs"),
                                        colon.hurst.rank %>% cbind(cell = "colon"),
                                        ileum.hurst.rank %>% cbind(cell = "ileum")) %>% filter(data %in% c("infected", "simulated")) %>%
  dplyr::select(-c(Hal, Ht, Hrs, Hs, `RSIboot FDR`)) %>% pivot_wider(names_from = data, values_from = He)

hurst.infected.vs.simulated.tb$group <- NA
hurst.infected.vs.simulated.tb[hurst.infected.vs.simulated.tb$infected < 0.5 | hurst.infected.vs.simulated.tb$simulated < 0.5, "group"] <- "anti-persistent"
hurst.infected.vs.simulated.tb[hurst.infected.vs.simulated.tb$infected >= 0.5 & hurst.infected.vs.simulated.tb$infected < 0.7 |
                                 hurst.infected.vs.simulated.tb$simulated >= 0.5 & hurst.infected.vs.simulated.tb$simulated < 0.7, "group"] <- "persistent"
hurst.infected.vs.simulated.tb[hurst.infected.vs.simulated.tb$infected >= 0.7 | hurst.infected.vs.simulated.tb$simulated >= 0.7, "group"] <- "Hurst phenomenon"
hurst.infected.vs.simulated.tb$group <- factor(hurst.infected.vs.simulated.tb$group, levels = c("anti-persistent", "persistent", "Hurst phenomenon"))

hurst.infected.vs.simulated.gg <- hurst.infected.vs.simulated.tb %>%
  ggplot(aes(x = infected, y = simulated, color = group)) + geom_point(size = 1) +
  scale_color_manual(values = c("grey40", "grey60", "tomato")) +
  facet_wrap(~factor(cell, level = c("hBECs", "colon", "ileum"))) +
  ylim(.25, 1) + xlim(.25, 1) +
  geom_density_2d(color = "black", linewidth = .25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = .25) +
  theme(legend.title = element_blank()) +
#  geom_vline(xintercept = .7, linewidth = .25, color = "red", linetype = "dashed") +
#  geom_hline(yintercept = .7, linewidth = .25, color = "red", linetype = "dashed") +
  labs(tag = "a")

# select genes above a He threshold

hurst.significant.genes <- rbind(hBECs.hurst.rank %>% cbind(cell = "hBECs"),
                                 colon.hurst.rank %>% cbind(cell = "colon"),
                                 ileum.hurst.rank %>% cbind(cell = "ileum")) %>% filter(data %in% c("infected", "simulated")) %>%
  dplyr::select(-c(Hal, Ht, Hrs, Hs, `RSIboot FDR`)) %>% pivot_wider(names_from = data, values_from = He) %>% 
  dplyr::filter(infected >= .7 & simulated < .7)

hurst.significant.genes.intersect <- hurst.significant.genes %>% dplyr::filter(gene %in%
                                                                                 names(table(hurst.significant.genes$gene)[table(hurst.significant.genes$gene) == 3]))

hurst.significant.genes.intersect.GO <- enrichGO(hurst.significant.genes.intersect$gene %>% unique(),
                                                 org.Hs.eg.db, keyType = "SYMBOL",
                                                 ont = "BP") %>% as.data.frame()

hurst.virus.process.targets <- hurst.significant.genes.intersect.GO["GO:0016032", "geneID"] %>% strsplit("/") %>% unlist()

# plot GO results

hurst.significant.genes.intersect.He <- lapply(setNames(nm = c("hBECs", "colon", "ileum")), function(x) {
  tb <- hurst.significant.genes.intersect %>% dplyr::filter(cell == x)
  lst <- tb$infected
  names(lst) <- tb$gene
  lst
})

hurst.significant.genes.intersect.bar.GO <- gg.bar.GO(hurst.significant.genes.intersect.GO[1:50,],
                                                      gene.values = hurst.significant.genes.intersect.He, legend = "median\nH") %>%
  annotate_figure(fig.lab = " b", fig.lab.size = 14)

# save plot

hurst.gg <- ggarrange(hurst.infected.vs.simulated.gg + theme(aspect.ratio = 1,
                                                             plot.tag.position = c(-.109, 1)),
                      hurst.significant.genes.intersect.bar.GO,
                      heights = c(1, 2.8),
                      ncol = 1)

ggsave("figs_out/final/hurst.eps", device = "eps",  hurst.gg, width = 6.85, height = 9)

# save Hurst and GO results

write_xlsx(list(hBECs_Hurst = hBECs.hurst.rank %>% dplyr::select(gene, data, He),
                colon_Hurst = colon.hurst.rank %>% dplyr::select(gene, data, He),
                ileum_Hurst = ileum.hurst.rank %>% dplyr::select(gene, data, He),

                GO_results = hurst.significant.genes.intersect.GO %>%
                  add.value.to.GO(hurst.significant.genes.intersect.He)), "Hurst_supplementary.xlsx")

# plot fit to power law

get.H <- function(hurst) hurst %>%
  dplyr::select(He) %>% unlist() %>% 
  setNames(nm = hurst %>% dplyr::select(gene) %>% unlist())

get.significant.H <- function(hurst) {
  
  significant.genes <- hurst %>% filter(data %in% c("infected", "simulated")) %>%
    dplyr::select(-c(Hal, Ht, Hrs, Hs)) %>% pivot_wider(names_from = data, values_from = He) %>% 
    filter(infected >= .7 & simulated < .7) %>% dplyr::select(gene)

  ifelse(unique(hurst$gene) %in% significant.genes$gene, TRUE, FALSE) %>% 
    setNames(nm = hurst %>% dplyr::select(gene) %>% unlist() %>% unique())
}

gg.fit.rank.hurst <- plot.power_law(list(hBECs = hBECs.TL$log_log.mean_sd$infected,
                                         colon = colon.TL$log_log.mean_sd$infected,
                                         ileum = ileum.TL$log_log.mean_sd$infected),
                                    list(hBECs = hBECs.hurst.rank %>% get.H(),
                                         colon = colon.hurst.rank %>% get.H(),
                                         ileum = ileum.hurst.rank %>% get.H()), "H", c(0, 100)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = .25) +
  geom_abline(intercept = -1.7, slope = .5, color = "black", linewidth = .25, linetype = "dashed") +
  scale_color_viridis() +
  xlab("mean (log)") + ylab("standard deviation (log)")

gg.fit.rank.hurst.significant <- plot.power_law(list(hBECs = hBECs.TL$log_log.mean_sd$infected,
                                                     colon = colon.TL$log_log.mean_sd$infected,
                                                     ileum = ileum.TL$log_log.mean_sd$infected),
                                                list(hBECs = hBECs.hurst.rank %>% get.significant.H(),
                                                     colon = colon.hurst.rank %>% get.significant.H(),
                                                     ileum = ileum.hurst.rank %>% get.significant.H()), "Ha", NULL) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = .25) +
  geom_abline(intercept = -1.7, slope = .5, color = "black", linewidth = .25, linetype = "dashed") +
  xlab("mean (log)") + ylab("standard deviation (log)") + theme(legend.position = "none") +
  scale_color_manual(values = c("grey60", "tomato"))

# mean vs H

hurst.mean.vs.hurst.tb <- hurst.infected.vs.simulated.tb
hurst.mean.vs.hurst.tb$`mean (log)` <- NA

hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "hBECs", "mean (log)"] <- hBECs.TL$log_log.mean_sd$infected$mean[match(hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "hBECs",]$gene,
                                                                                                                             rownames(hBECs.TL$log_log.mean_sd$infected))] %>% exp() %>% log10()
hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "colon", "mean (log)"] <- colon.TL$log_log.mean_sd$infected$mean[match(hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "colon",]$gene,
                                                                                                                             rownames(colon.TL$log_log.mean_sd$infected))] %>% exp() %>% log10()
hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "ileum", "mean (log)"] <- ileum.TL$log_log.mean_sd$infected$mean[match(hurst.mean.vs.hurst.tb[hurst.mean.vs.hurst.tb$cell == "ileum",]$gene,
                                                                                                                             rownames(ileum.TL$log_log.mean_sd$infected))] %>% exp() %>% log10()

hurst.mean.vs.h.gg <- hurst.mean.vs.hurst.tb %>%
  ggplot(aes(`mean (log)`, infected)) + geom_point() + facet_wrap(~cell) +
  geom_hline(yintercept = .7, linewidth = .25, color = "red", linetype = "dashed") +
  geom_hline(yintercept = .5, linewidth = .25, color = "red", linetype = "solid")

# save H values + power law figs

rank.hurst.taylor.gg <- ggarrange(hurst.all_data.gg + labs(tag = "a") + ylab("distribution"),
                                  gg.fit.rank.hurst + labs(tag = "b"),
                                  hurst.mean.vs.h.gg + labs(tag = "c") + ylab("H"),
                                  ncol = 1)

ggsave("figs_out/final/rank_hurst_taylor.eps", device = "eps",  rank.hurst.taylor.gg, width = 6.85, height = 7.5)

# test if H is higher in any cell for all genes

unpaired.wilcox.hurst.all.genes <- list(`hBECs vs colon` = wilcox.test(hBECs.hurst.rank[hBECs.hurst.rank$data == "infected",]$He,
                                                                       colon.hurst.rank[colon.hurst.rank$data == "infected",]$He,
                                                                       paired = FALSE),
                                        `hBECs vs ileum` = wilcox.test(hBECs.hurst.rank[hBECs.hurst.rank$data == "infected",]$He,
                                                                       ileum.hurst.rank[ileum.hurst.rank$data == "infected",]$He,
                                                                       paired = FALSE),
                                        `ileum vs colon` = wilcox.test(ileum.hurst.rank[ileum.hurst.rank$data == "infected",]$He,
                                                                       colon.hurst.rank[colon.hurst.rank$data == "infected",]$He,
                                                                       paired = FALSE))

unpaired.wilcox.hurst.all.genes.p.adj <- unpaired.wilcox.hurst.all.genes %>% lapply(function(x) x$p.value) %>% unlist() %>% p.adjust(method = "fdr")

# test if H is higher in any cell for significant H intersect genes

paired.wilcox.hurst.significant.genes.intersect <- list(`hBECs vs colon` = wilcox.test(hurst.significant.genes.intersect.He$hBECs[sort(names(hurst.significant.genes.intersect.He$hBECs))],
                                                                                       hurst.significant.genes.intersect.He$colon[sort(names(hurst.significant.genes.intersect.He$colon))],
                                                                                       paired = TRUE),
                                                        `hBECs vs ileum` = wilcox.test(hurst.significant.genes.intersect.He$hBECs[sort(names(hurst.significant.genes.intersect.He$hBECs))],
                                                                                       hurst.significant.genes.intersect.He$ileum[sort(names(hurst.significant.genes.intersect.He$ileum))],
                                                                                       paired = TRUE),
                                                        `ileum vs colon` = wilcox.test(hurst.significant.genes.intersect.He$ileum[sort(names(hurst.significant.genes.intersect.He$ileum))],
                                                                                       hurst.significant.genes.intersect.He$colon[sort(names(hurst.significant.genes.intersect.He$colon))],
                                                                                       paired = TRUE))

paired.wilcox.hurst.significant.genes.intersect.p.adj <- paired.wilcox.hurst.significant.genes.intersect %>% lapply(function(x) x$p.value) %>% unlist() %>% p.adjust(method = "fdr")

# data structure figure, will be latter edited

struc.theme <- ttheme_default(base_size = 7)

hBECs.abundance.small <- hBECs.mtx$matrices$infected %>% normalize.cells() %>% as.data.frame()
hBECs.abundance.small <- hBECs.abundance.small[c("EIF3A", "RPL27", "A1BG"), 1:5] %>% round(5)
rownames(hBECs.abundance.small) <- c("gene a", "gene b", "gene c")
colnames(hBECs.abundance.small) <- c("cell i", "cell i + 1", "cell 100", "cell 101", "cell 200")
hBECs.abundance.small$... <- ""
hBECs.abundance.small["...",] <-  ""

struc.a <- ggarrange(grid.arrange(tableGrob(hBECs.abundance.small[c(1, 2, 6)], theme = struc.theme)),
                     ggplot(hBECs.TL$log_log.mean_sd$infected %>% mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)),
                            aes(mean, sd)) + theme(aspect.ratio = 1) + ggtitle("") + # just to scale better
                       geom_point(color = "grey60") + xlab("mean (log)") + ylab("standard deviation (log)")) %>% annotate_figure(fig.lab = "a", fig.lab.size = 14)

colnames(hBECs.abundance.small)[1] <- "cell 1"

struc.b <- ggarrange(grid.arrange(tableGrob(hBECs.abundance.small[c(1, 6, 3)], theme = struc.theme) %>% grid.arrange(top = text_grob("bin i", size = 13, vjust = 1.5)),
                                  tableGrob(hBECs.abundance.small[c(4, 6, 5)], theme = struc.theme) %>% grid.arrange(top = text_grob("bin i + 1", size = 13, vjust = 1.5)), #bottom = "viral load bin i < viral load bin i + 1",
                                  ncol = 1),
                     ggarrange(ggplot(hBECs.bin30.TL$log_log.mean_sd[[1]] %>% mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)),
                                         aes(mean, sd)) +
                                 theme(aspect.ratio = 1) + ggtitle("bin i") +
                                 geom_point(color = "grey60") + xlab("mean (log)") + ylab("standard deviation (log)"),
                                 ggplot(hBECs.bin30.TL$log_log.mean_sd[[2]] %>% mutate(mean = log(exp(mean), 10), sd = log(exp(sd), 10)),
                                        aes(mean, sd)) +
                                   theme(aspect.ratio = 1) + ggtitle("bin i + 1") +
                                 geom_point(color = "grey60") + xlab("mean (log)") + ylab("standard deviation (log)"),
                               ncol = 1)) %>% annotate_figure(fig.lab = "b", fig.lab.size = 14)

hBECs.mean.small <- hBECs.bin30_mean_sd.mtx$mean[c("EIF3A", "RPL27", "A1BG"), 1:2] %>% round(5)
rownames(hBECs.mean.small) <- c("gene a", "gene b", "gene c")
colnames(hBECs.mean.small) <- c("cells 1-100", "cells 101-200")
hBECs.mean.small$... <- ""
hBECs.mean.small["...",] <-  ""

hBECs.rank.small <- hBECs.rank$infected %>% as.data.frame()
hBECs.rank.small <- hBECs.rank.small[c("EIF3A", "RPL27", "A1BG"), 1:4]
rownames(hBECs.rank.small) <- c("gene a", "gene b", "gene c")
colnames(hBECs.rank.small) <- c("cell i", "cell i + 1", "cells 1-100", "cells 101-200")
hBECs.rank.small$... <-  ""
hBECs.rank.small["...",] <-  ""

struc.c.d <- ggarrange(grid.arrange(tableGrob(hBECs.rank.small[c(1, 2, 5)], theme = struc.theme), top = "") %>%
                                      annotate_figure(fig.lab = "c", fig.lab.size = 14),
                       ggarrange(grid.arrange(tableGrob(hBECs.mean.small[c(1, 2, 3)], theme = struc.theme),
                                              top = text_grob("mean gene abundances", size = 10, vjust = 1)),
                                 grid.arrange(tableGrob(hBECs.rank.small[3:5], theme = struc.theme),
                                              top = text_grob("rank of mean gene abundances", size = 10, vjust = 1))) %>%
                         annotate_figure(fig.lab = "d", fig.lab.size = 14),
                       widths = c(1, 2.5), nrow = 1)

struc.gg <- ggarrange(struc.a, struc.b, struc.c.d, ncol = 1, heights = c(1, 2, 1))

ggsave("figs_out/final/data_structure.eps", struc.gg, device = "eps",  height = 8, width = 6.85)

# network analyses

huri_edges <- read.table("HuRI.tsv")
huri <- graph_from_edgelist(huri_edges %>% as.matrix(), directed = FALSE)

huri.k <- degree(huri)
huri.b <- betweenness(huri, directed = FALSE)

# hubs based on 1.5xIQR rule and bottlenecks = articulation points

huri.k.1.5_iqr <- IQR(huri.k) * 1.5
huri.k.1.5_iqr.hubs <- huri.k[huri.k > huri.k.1.5_iqr]

# convert symbol to ensembl

RSI.boot.bin30_intersect.ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = RSI.boot.bin30_intersect, columns = "ENSEMBL", keytype = "SYMBOL")

hurst.significant.genes.intersect.genes <- hurst.significant.genes.intersect$gene %>% unique()
hurst.significant.genes.intersect.ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = hurst.significant.genes.intersect.genes, columns = "ENSEMBL", keytype = "SYMBOL")

# hubs

RSI.boot.bin30_intersect.hubs <- huri.k.1.5_iqr.hubs[RSI.boot.bin30_intersect.ensembl$ENSEMBL[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(huri.k.1.5_iqr.hubs)]]
RSI.boot.bin30_intersect.hubs.symbol <- RSI.boot.bin30_intersect.ensembl[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(RSI.boot.bin30_intersect.hubs),]$SYMBOL

hurst.significant.genes.intersect.hubs <- huri.k.1.5_iqr.hubs[hurst.significant.genes.intersect.ensembl$ENSEMBL[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(huri.k.1.5_iqr.hubs)]]
hurst.significant.genes.intersect.hubs.symbol <- hurst.significant.genes.intersect.ensembl[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(hurst.significant.genes.intersect.hubs),]$SYMBOL

# articulation points

huri.articulation_points <- articulation_points(huri) %>% as.list() %>% unlist()

RSI.boot.bin30_intersect.ap <- huri.articulation_points[RSI.boot.bin30_intersect.ensembl$ENSEMBL[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(huri.articulation_points)]]
RSI.boot.bin30_intersect.ap.symbol <- RSI.boot.bin30_intersect.ensembl[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(RSI.boot.bin30_intersect.ap),]$SYMBOL

hurst.significant.genes.intersect.ap <- huri.articulation_points[hurst.significant.genes.intersect.ensembl$ENSEMBL[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(huri.articulation_points)]]
hurst.significant.genes.intersect.ap.symbol <- hurst.significant.genes.intersect.ensembl[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(hurst.significant.genes.intersect.ap),]$SYMBOL

# save network results

write_xlsx(list(hubs = AnnotationDbi::select(org.Hs.eg.db, keys = names(huri.k.1.5_iqr.hubs), columns = "SYMBOL", keytype = "ENSEMBL"),
                `bottlenecks (ap)` = AnnotationDbi::select(org.Hs.eg.db, keys = names(huri.articulation_points), columns = "SYMBOL", keytype = "ENSEMBL"),
                
                `common punctual stab hubs` = RSI.boot.bin30_intersect.ensembl[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(RSI.boot.bin30_intersect.hubs),],
                `common punctual stab ap` = RSI.boot.bin30_intersect.ensembl[RSI.boot.bin30_intersect.ensembl$ENSEMBL %in% names(RSI.boot.bin30_intersect.ap),],
                
                `common persistent rank hubs` = hurst.significant.genes.intersect.ensembl[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(hurst.significant.genes.intersect.hubs),],
                `common persistent rank ap` = hurst.significant.genes.intersect.ensembl[hurst.significant.genes.intersect.ensembl$ENSEMBL %in% names(hurst.significant.genes.intersect.ap),]),
           "Network_supplementary.xlsx")

# genes not detected previously in the used datasets

hurst.RSI.intersect.significant <- c(hurst.significant.genes.intersect.genes, RSI.boot.bin30_intersect) %>% unique()

hBECs.ori <- read_xlsx("prev_res/journal.pbio.3001143.s006.xlsx")
hBECs.ori.DE <- hBECs.ori %>% dplyr::filter(pval_corrected < .05)
hBECs.ori.not_DE <- hBECs.ori %>% dplyr::filter(pval_corrected >= .05)

intersect_genes_not_in_hBECs.ori <- hBECs.ori.not_DE[hBECs.ori.not_DE$Gene %in% hurst.RSI.intersect.significant,]$Gene %>% unique() %>% sort()

load("prev_res/COVID19_July.rda")

# do tests 24h vs mock and 24h vs bystander with MAST

Idents(Colon_H_T) <- Colon_H_T$Groups
Colon_H_T.DE.ALL <- FindMarkers(Colon_H_T, ident.1 = "24h_Infected", ident.2 = "Mock_Non-Infected", test.use = "MAST", assay = "RNA")

Colon_H_T$CellTypes.Groups <- paste0(Colon_H_T$CellTypes, "_", Colon_H_T$Groups)
Idents(Colon_H_T) <- Colon_H_T$CellTypes.Groups
Colon_H_T.DE.IE2 <- FindMarkers(Colon_H_T, ident.1 = "Inmature Enterocyte 2_24h_Infected", ident.2 = "Inmature Enterocyte 2_Mock_Non-Infected", test.use = "MAST", assay = "RNA")
Colon_H_T.DE.IE2.bystander <- FindMarkers(Colon_H_T, ident.1 = "Inmature Enterocyte 2_24h_Bystander", ident.2 = "Inmature Enterocyte 2_Mock_Non-Infected", test.use = "MAST", assay = "RNA")
Colon_H_T.DE.IE2.vs_bystander <- FindMarkers(Colon_H_T, ident.1 = "Inmature Enterocyte 2_24h_Infected", ident.2 = "Inmature Enterocyte 2_24h_Bystander", test.use = "MAST", assay = "RNA")

Idents(Illeum_H_T) <- Illeum_H_T$Groups
Illeum_H_T.DE.ALL <- FindMarkers(Illeum_H_T, ident.1 = "24h_Infected", ident.2 = "Mock_Non-Infected", test.use = "MAST", assay = "RNA")

Illeum_H_T$CellTypes.Groups <- paste0(Illeum_H_T$CellTypes, "_", Illeum_H_T$Groups)
Idents(Illeum_H_T) <- Illeum_H_T$CellTypes.Groups
Illeum_H_T.DE.IE2 <- FindMarkers(Illeum_H_T, ident.1 = "Inmature Enterocyte 2_24h_Infected", ident.2 = "Inmature Enterocyte 2_Mock_Non-Infected", test.use = "MAST", assay = "RNA")
Illeum_H_T.DE.IE2.bystander <- FindMarkers(Illeum_H_T, ident.1 = "Inmature Enterocyte 2_24h_Bystander", ident.2 = "Inmature Enterocyte 2_Mock_Non-Infected", test.use = "MAST", assay = "RNA")
Illeum_H_T.DE.IE2.vs_bystander <- FindMarkers(Illeum_H_T, ident.1 = "Inmature Enterocyte 2_24h_Infected", ident.2 = "Inmature Enterocyte 2_24h_Bystander", test.use = "MAST", assay = "RNA")

ALL.hIECs.DE <- c(rownames(Colon_H_T.DE.ALL), rownames(Colon_H_T.DE.IE2), rownames(Colon_H_T.DE.IE2.bystander), rownames(Colon_H_T.DE.IE2.vs_bystander),
                  rownames(Illeum_H_T.DE.ALL), rownames(Illeum_H_T.DE.IE2), rownames(Illeum_H_T.DE.IE2.bystander), rownames(Illeum_H_T.DE.IE2.vs_bystander)) %>% unique()

intersect_genes_not_in_hIECs.ori <- hurst.RSI.intersect.significant[!hurst.RSI.intersect.significant %in% ALL.hIECs.DE]

significant.genes.this.study <- intersect_genes_not_in_hBECs.ori[intersect_genes_not_in_hBECs.ori %in% intersect_genes_not_in_hIECs.ori]
significant.genes.this.study.ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = significant.genes.this.study, columns = "ENSEMBL", keytype = "SYMBOL")

significant.genes.this.study.hubs <- significant.genes.this.study.ensembl[significant.genes.this.study.ensembl$ENSEMBL %in% names(huri.k.1.5_iqr.hubs),]$SYMBOL
significant.genes.this.study.ap <- significant.genes.this.study.ensembl[significant.genes.this.study.ensembl$ENSEMBL %in% names(huri.articulation_points),]$SYMBOL

# save genes detected in this study and not on the original ones

write_xlsx(list(`significant genes all cells` = significant.genes.this.study,
                `hubs` = significant.genes.this.study.hubs,
                `bottlenecks` = significant.genes.this.study.ap) %>% lapply(as.data.frame),
           "Significant_genes_this_study_supplementary.xlsx", col_names = FALSE)

# ICD based on 10.1080/2162402X.2015.1069938

ICD.genes <- c("ENTPD1", "NT5E", "CALR", "HMGB1", "HSP90AA1",
               "ATG5", "BAX", "CASP8", "PDIA3", "EIF2AK3", "PIK3CA",
               "CXCR3", "IFNA1", "IFNB1", "IL10", "IL6", "TNF",
               "CASP1", "IL1R1", "IL1B", "NLRP3", "P2RX7",
               "LY96", "MYD88", "TLR4",
               "CD4", "CD8A", "CD8B", "FOXP3",
               "IFNG", "IFNGR1", "IL17A", "IL17RA", "PRF1")

hurst.RSI.intersect.significant[hurst.RSI.intersect.significant %in% ICD.genes]
RSI.boot.bin30_intersect[RSI.boot.bin30_intersect %in% ICD.genes]

# apopotisis-related genes in the interactome

RSI.boot.bin30_intersect.apoptosis <- RSI.boot.bin30_intersect.GO[grepl("apopto|p53", RSI.boot.bin30_intersect.GO$Description),]$geneID %>%
  lapply(strsplit, "/") %>% unlist() %>% unique()

hurst.significant.genes.intersect.apoptosis <- hurst.significant.genes.intersect.GO[grepl("apopto|p53", hurst.significant.genes.intersect.GO$Description),]$geneID %>%
  lapply(strsplit, "/") %>% unlist() %>% unique()

RSI.boot.bin30_intersect.apoptosis.ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = RSI.boot.bin30_intersect.apoptosis, columns = "ENSEMBL", keytype = "SYMBOL")
RSI.boot.bin30_intersect.apoptosis.hubs <- RSI.boot.bin30_intersect.apoptosis.ensembl[RSI.boot.bin30_intersect.apoptosis.ensembl$ENSEMBL %in% names(huri.k.1.5_iqr.hubs),]$SYMBOL
RSI.boot.bin30_intersect.apoptosis.ap <- RSI.boot.bin30_intersect.apoptosis.ensembl[RSI.boot.bin30_intersect.apoptosis.ensembl$ENSEMBL %in% names(huri.articulation_points),]$SYMBOL

hurst.significant.genes.intersect.apoptosis.ensemble <-AnnotationDbi::select(org.Hs.eg.db, keys = hurst.significant.genes.intersect.apoptosis, columns = "ENSEMBL", keytype = "SYMBOL")
hurst.significant.genes.intersect.apoptosis.hubs <- hurst.significant.genes.intersect.apoptosis.ensemble[hurst.significant.genes.intersect.apoptosis.ensemble$ENSEMBL %in% names(huri.k.1.5_iqr.hubs),]$SYMBOL
hurst.significant.genes.intersect.apoptosis.ap <- hurst.significant.genes.intersect.apoptosis.ensemble[hurst.significant.genes.intersect.apoptosis.ensemble$ENSEMBL %in% names(huri.articulation_points),]$SYMBOL

# interactions of bottlenecks

RSI.boot.bin30_intersect.ap.interactions.GO <- setNames(nm = RSI.boot.bin30_intersect.ap.symbol) %>%
  lapply(function(x) {
    ens <- AnnotationDbi::select(org.Hs.eg.db, keys = x, columns = "ENSEMBL", keytype = "SYMBOL")$ENSEMBL
    interactions.ens <- huri_edges[huri_edges$V1 %in% ens | huri_edges$V2 %in% ens,] %>% unlist() %>% unique()
    interactions.ens <- interactions.ens[!interactions.ens %in% ens]
    interactions.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = interactions.ens, columns = "SYMBOL", keytype = "ENSEMBL")$SYMBOL
    GO <- enrichGO(interactions.symbol,
                   org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP") %>% as.data.frame()
    list(interactions = interactions.symbol,
         GO = GO)
  })

hurst.significant.genes.intersect.ap.interactions.GO <- setNames(nm = hurst.significant.genes.intersect.ap.symbol) %>%
  lapply(function(x) {
    ens <- AnnotationDbi::select(org.Hs.eg.db, keys = x, columns = "ENSEMBL", keytype = "SYMBOL")$ENSEMBL
    interactions.ens <- huri_edges[huri_edges$V1 %in% ens | huri_edges$V2 %in% ens,] %>% unlist() %>% unique()
    interactions.ens <- interactions.ens[!interactions.ens %in% ens]
    interactions.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = interactions.ens, columns = "SYMBOL", keytype = "ENSEMBL")$SYMBOL
    GO <- enrichGO(interactions.symbol,
                   org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP") %>% as.data.frame()
    list(interactions = interactions.symbol,
         GO = GO)
  })

