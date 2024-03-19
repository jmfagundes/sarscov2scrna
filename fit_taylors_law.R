source("source_me.R")

# prepare gene matrices

hBECs.mtx <- process.gene.count("count_matrix_hBECs.tbl", remove.MT = remove.MT,
                                min.rna = list(`0dpi` = 1878, `1dpi` = 1232, `2dpi` = 1752, `3dpi` = 1781),
                                max.rna = list(`0dpi` = 16810, `1dpi` = 30555, `2dpi` = 20120, `3dpi` = 13181),
                                raw.mtx = TRUE)
colon.mtx <- process.gene.count("count_matrix_colon.tbl", remove.MT = remove.MT,
                                min.rna = list(`0hpi` = 1030, `12hpi` = 1149, `24hpi` = 1232),
                                max.rna = list(`0hpi` = 47489, `12hpi` = 46952, `24hpi` = 31513),
                                raw.mtx = TRUE)
ileum.mtx <- process.gene.count("count_matrix_ileum.tbl", remove.MT = remove.MT,
                                min.rna = list(`0hpi` = 968, `12hpi` = 1100, `24hpi` = 908),
                                max.rna = list(`0hpi` = 40234, `12hpi` = 40163, `24hpi` = 42604),
                                raw.mtx = TRUE)

# calculate mean viral counts in empty droplets

hBECs.empty_droplets <- list(`0dpi` = hBECs.mtx$raw.matrices$`0dpi`[colSums(hBECs.mtx$raw.matrices$`0dpi`) < 1878],
                             `1dpi` = hBECs.mtx$raw.matrices$`1dpi`[colSums(hBECs.mtx$raw.matrices$`1dpi`) < 1232],
                             `2dpi` = hBECs.mtx$raw.matrices$`2dpi`[colSums(hBECs.mtx$raw.matrices$`2dpi`) < 1752],
                             `3dpi` = hBECs.mtx$raw.matrices$`3dpi`[colSums(hBECs.mtx$raw.matrices$`3dpi`) < 1781])

colon.empty_droplets <- list(`0hpi` = colon.mtx$raw.matrices$`0hpi`[colSums(colon.mtx$raw.matrices$`0hpi`) < 1030],
                             `12hpi` = colon.mtx$raw.matrices$`12hpi`[colSums(colon.mtx$raw.matrices$`12hpi`) < 1149],
                             `24hpi` = colon.mtx$raw.matrices$`24hpi`[colSums(colon.mtx$raw.matrices$`24hpi`) < 1232])

ileum.empty_droplets <- list(`0hpi` = ileum.mtx$raw.matrices$`0hpi`[colSums(ileum.mtx$raw.matrices$`0hpi`) < 968],
                             `12hpi` = ileum.mtx$raw.matrices$`12hpi`[colSums(ileum.mtx$raw.matrices$`12hpi`) < 1100],
                             `24hpi` = ileum.mtx$raw.matrices$`24hpi`[colSums(ileum.mtx$raw.matrices$`24hpi`) < 908])

# identify infected cells (> mean ambient viral RNA + 2 SD)
# if no viral RNA in empty droplets, set threshold to 10 (which was always the case)

hBECs.V.min <- hBECs.empty_droplets %>% lapply(function(x) {
  threshold <- mean(x["SARS_Cov2",]) + (sd(x["SARS_Cov2",]) * 2)
  if (is.na(threshold)) return(10)
  else return(threshold)
})

colon.V.min <- colon.empty_droplets %>% lapply(function(x) {
  threshold <- mean(x["SARS_Cov2",]) + (sd(x["SARS_Cov2",]) * 2)
  if (is.na(threshold)) return(10)
  else return(threshold)
})

ileum.V.min <- ileum.empty_droplets %>% lapply(function(x) {
  threshold <- mean(x["SARS_Cov2",]) + (sd(x["SARS_Cov2",]) * 2)
  if (is.na(threshold)) return(10)
  else return(threshold)
})

# load TPT data

hBECs.mtx.TPT <- process.gene.count("Expected_TPT_matrix_hBECs_final.csv", remove.MT = remove.MT, TPT = TRUE,
                                    infected.cells = names(hBECs.mtx$viral_accumulation$infected))
colon.mtx.TPT <- process.gene.count("Expected_TPT_matrix_colon_final.csv", remove.MT = remove.MT, TPT = TRUE,
                                    infected.cells = names(colon.mtx$viral_accumulation$infected))
ileum.mtx.TPT <- process.gene.count("Expected_TPT_matrix_ileum_final.csv", remove.MT = remove.MT, TPT = TRUE,
                                    infected.cells = names(ileum.mtx$viral_accumulation$infected))

# fit to Taylor's power law

hBECs.TL <- fit.TL(hBECs.mtx$matrices, prefix.xlsx = NULL)
colon.TL <- fit.TL(colon.mtx$matrices, prefix.xlsx = NULL)
ileum.TL <- fit.TL(ileum.mtx$matrices, prefix.xlsx = NULL)

hBECs.TL.TPT <- fit.TL(hBECs.mtx.TPT$matrices, prefix.xlsx = NULL)
colon.TL.TPT <- fit.TL(colon.mtx.TPT$matrices, prefix.xlsx = NULL)
ileum.TL.TPT <- fit.TL(ileum.mtx.TPT$matrices, prefix.xlsx = NULL)

# Taylor's parameters

taylors.params <- rbind(cbind(hBECs.TL$params, cell = "hBECs"),
                        cbind(colon.TL$params, cell = "colon"),
                        cbind(ileum.TL$params, cell = "ileum"))

taylors.params.TPT <- rbind(cbind(hBECs.TL.TPT$params, cell = "hBECs"),
                            cbind(colon.TL.TPT$params, cell = "colon"),
                            cbind(ileum.TL.TPT$params, cell = "ileum"))

taylors.params$cell <- factor(taylors.params$cell, levels = c("hBECs", "colon", "ileum"))
taylors.params$data <- taylors.params$data %>% gsub("uninfected_0dpi", "uninfected 0dpi", .) %>% gsub("uninfected_", "bystander ", .) %>% gsub("_", " ", .)
colnames(taylors.params) <- c("data", "V", "beta", "R^2", "number of cells", "number of genes", "sparsity", "model", "cell")

taylors.params.TPT$cell <- factor(taylors.params.TPT$cell, levels = c("hBECs", "colon", "ileum"))
taylors.params.TPT$data <- taylors.params.TPT$data %>% gsub("uninfected_0dpi", "uninfected 0dpi", .) %>% gsub("uninfected_", "bystander ", .) %>% gsub("_", " ", .)
colnames(taylors.params.TPT) <- c("data", "V", "beta", "R^2", "number of cells", "number of genes", "sparsity", "model", "cell")

gg.params <- ggplot(taylors.params, aes(x = V, y = beta, label = data, color = `R^2`, size = `number of cells`)) + geom_point() +
  facet_wrap(~cell) +
  geom_text_repel(seed = 12345, size = 1.22, segment.size = 0.25, force = 75, max.overlaps = Inf) +
  theme(legend.key.size = unit(.3, 'cm'), legend.text = element_text(size = 6), legend.title = element_text(size = 8)) +
  scale_radius(range = c(.5, 2))

gg.params.TPT <- ggplot(taylors.params.TPT, aes(x = V, y = beta, label = data, color = `R^2`, size = `number of cells`)) + geom_point() +
  facet_wrap(~cell) +
  geom_text_repel(seed = 12345, size = 1.22, segment.size = 0.25, force = 75, max.overlaps = Inf) +
  theme(legend.key.size = unit(.3, 'cm'), legend.text = element_text(size = 6), legend.title = element_text(size = 8)) +
  scale_radius(range = c(.5, 2))

# save R objects for other scripts

dir.create("RData")
dir.create(paste0("RData/MT"))
dir.create(paste0("RData/TPT_MT"))

save(hBECs.mtx, colon.mtx, ileum.mtx,
     hBECs.TL, colon.TL, ileum.TL,
     file = "RData/MT/cell_infection_pseudotime_matrices_TL.RData")

save(hBECs.mtx.TPT, colon.mtx.TPT, ileum.mtx.TPT,
     hBECs.TL.TPT, colon.TL.TPT, ileum.TL.TPT,
     file = "RData/TPT_MT/cell_infection_pseudotime_matrices_TL.RData")

# plot plots

dir.create("figs_out")
dir.create("figs_out/MT")
dir.create("figs_out/TPT_MT")

ggsave("figs_out/MT/taylors_parameters.pdf", gg.params, width = 6.85, height = 2.5)
ggsave("figs_out/TPT_MT/taylors_parameters.pdf", gg.params.TPT, width = 6.85, height = 2.5)

# simplified figs

gg.params.TPT <- ggplot(taylors.params.TPT %>% dplyr::filter(data == "infected"), aes(x = V, y = beta, color = cell)) + geom_point(size = 5) +
  theme(legend.key.size = unit(.3, 'cm'), legend.text = element_text(size = 6), legend.title = element_text(size = 8)) +
  scale_radius(range = c(.5, 2))

ggsave("figs_out/TPT_MT/taylors_parameters_simplified.png", gg.params.TPT, width = 6.85, height = 6, dpi = 300)

# segmented fit

hBECs.seg.TL <- fit.seg.TL(hBECs.TL)
colon.seg.TL <- fit.seg.TL(colon.TL)
ileum.seg.TL <- fit.seg.TL(ileum.TL)

hBECs.seg.TL.TPT <- fit.seg.TL(hBECs.TL.TPT)
colon.seg.TL.TPT <- fit.seg.TL(colon.TL.TPT)
ileum.seg.TL.TPT <- fit.seg.TL(ileum.TL.TPT)

# F-test of the residuals

f.test.seg.vs.unseg <- list(hBECs = var.test(hBECs.seg.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             hBECs.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             alternative = "less"),
                            colon = var.test(colon.seg.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             colon.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             alternative = "less"),
                            ileum = var.test(ileum.seg.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             ileum.TL.TPT$lm.log_log.mean_sd$infected$residuals,
                                             alternative = "less"))

# infected cells parameters

taylors.params.TPT.infected <- taylors.params.TPT %>% dplyr::filter(data == "infected")

seg.fit.taylor.params.params.infected <- rbind(cbind(hBECs.seg.TL.TPT$params, cell = "hBECs"),
                                               cbind(colon.seg.TL.TPT$params, cell = "colon"),
                                               cbind(ileum.seg.TL.TPT$params, cell = "ileum")) %>% dplyr::filter(data == "infected")
