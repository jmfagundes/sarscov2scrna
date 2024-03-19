library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(readxl)
library(segmented)
library(ggpubr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(viridis)
library(ggrepel)
library(pracma)
library(grid)
library(gridExtra)

# packages that are not used in the final analysis, but some functions may use them for other types of data or analysis
#library(minpack.lm)
#library(statcomp)
#library(randtests)

# variables to control analysis, library of functions and other useful things

# variables

remove.MT <- FALSE
if (remove.MT) out.dir <- "noMT/" else out.dir <- "MT/"

use.TPT <- TRUE
if (use.TPT) out.dir <- paste0("TPT_", out.dir)

use.seg.fit <- TRUE # use segmented fit to get residuals

# control number of bins

hBECs.n.bins <- 30
colon.n.bins <- 30
ileum.n.bins <- 30

# functions

# return gene abundances

#' normalize.cells
#' 
#' Calculate gene abundances
#' 
#' @param mtx A matrix where rows are genes and columns are cells
#' 
#' @return output Gene abundances matrix
#' @export
#' 
#' @examples
normalize.cells <- function(mtx) t(t(mtx)/colSums(mtx))

#' process.gene.count
#' 
#' Process gene count matrices from Gutiérrez et al.
#' 
#' @param path Path to .tbl file
#' @param V.ann Virus annotation on the gene matrix
#' @param V.min Integer or named list of minimum virus counts per dataset
#' @param pct.infected.cells Get only the infected cells with the highest viral load to reach this infection rate, not used
#' @param min.cells Minimum number of cells to return a matrix
#' @param remove.MT Remove mitochondiral genes
#' @param infected.cells Mark these cells as infected
#' @param min.rna Vector of minimum number of transcripts in a cell per dataset
#' @param max.rna Vector of maximum number of transcripts in a cells per dataset
#' @param raw.mtx Return raw matrices
#' @param TPT Mark true if input is TPT
#' 
#' @return list List of a list of the output matrices and viral accumulation per cell per dataset
#' @export
#' 
#' @examples
process.gene.count <- function(path,
                               V.ann = "SARS_Cov2",
                               V.min = 10,
                               pct.infected.cells = NULL,
                               min.cells = 50,
                               remove.MT = FALSE,
                               infected.cells = NULL,
                               min.rna = NULL,
                               max.rna = NULL,
                               raw.mtx = FALSE,
                               TPT = FALSE) {
  
  matrices <- list()
  raw.matrices <- list()
  viral_accumulation <- list()
  
  if (!TPT) {
    gene.mtx <- read.table(path, header = TRUE, check.names = FALSE)
    rownames(gene.mtx) <- gene.mtx$cell_id
    gene.mtx <- gene.mtx[-1]
  } else {
    gene.mtx <- read.table(path, header = TRUE, check.names = FALSE, sep = ",", row.names = 1)
    gene.mtx <- gene.mtx[!colnames(gene.mtx) %in% c("total_counts", "log2_SARS", "log2_Astrovirus")]
  }
  
  # get datasets names
  
  datasets <- rownames(gene.mtx) %>% gsub("_.*", "", .) %>% unique() 
  
  # process each dataset
  
  for (d in datasets) {
    
    gene.mtx.dataset <- gene.mtx[grep(d, rownames(gene.mtx)),]
    if (raw.mtx) raw.matrices[[d]] <- gene.mtx.dataset %>% t() %>% as.data.frame()
    
    # get cell-containing GEMs
    
    if (!is.null(min.rna)) gene.mtx.dataset <- gene.mtx.dataset[rowSums(gene.mtx.dataset) >= min.rna[[d]],]
    if (!is.null(max.rna)) gene.mtx.dataset <- gene.mtx.dataset[rowSums(gene.mtx.dataset) <= max.rna[[d]],]
    
    n.cells <- gene.mtx[grep(d, rownames(gene.mtx)),] %>% nrow()
    
    message(paste0("raw number of cells in ", d, ": ", n.cells))
    message(paste0("number of cell-containing droplets: ", nrow(gene.mtx.dataset)))
    
    # get infected cells
    
    if (!is.null(infected.cells)) gene.mtx.infected <- gene.mtx.dataset[rownames(gene.mtx.dataset) %in% infected.cells,]
    else if (is.list(V.min)) gene.mtx.infected <- gene.mtx.dataset[gene.mtx.dataset[[V.ann]] >= V.min[[d]],]
    else gene.mtx.infected <- gene.mtx.dataset[gene.mtx.dataset[[V.ann]] >= V.min,]
    
    n.cells.vir.RNA <- nrow(gene.mtx.dataset[gene.mtx.dataset[[V.ann]] > 0,])
    
    # reorder infected cells and select uninfected
    # also transpose matrix so that 
    
    gene.mtx.infected <- gene.mtx.infected[order(gene.mtx.infected[[V.ann]]),] %>% t() 
    gene.mtx.uninfected <- gene.mtx.dataset[gene.mtx.dataset[[V.ann]] == 0,] %>% t()
    
    # get viral accumulation
    
    viral_accumulation.dataset <- gene.mtx.infected[V.ann,] / colSums(gene.mtx.infected)

    # remove mitochondrial genes and viral counts
    
    if (remove.MT) {
      MT.genes <- colnames(gene.mtx)[grep("^MT-", colnames(gene.mtx))]
      gene.mtx.infected <- gene.mtx.infected[!rownames(gene.mtx.infected) %in% c(V.ann, MT.genes),]
      gene.mtx.uninfected <- gene.mtx.uninfected[!rownames(gene.mtx.uninfected) %in% c(V.ann, MT.genes),]
    } else {
      gene.mtx.infected <- gene.mtx.infected[rownames(gene.mtx.infected) != V.ann,]
      gene.mtx.uninfected <- gene.mtx.uninfected[rownames(gene.mtx.uninfected) != V.ann,]
    }
    
    gene.mtx.infected <- gene.mtx.infected %>% as.data.frame()
    gene.mtx.uninfected <- gene.mtx.uninfected %>% as.data.frame()
    
    message(paste0("raw proportion of infected cells: ", ncol(gene.mtx.infected) / nrow(gene.mtx.dataset)))
    
    # filter infected cells based on expected infection rate, used only for some exploration of the data
    
    if (!is.null(pct.infected.cells)) {
      
      n.infected.cells <- n.cells * pct.infected.cells
      viral_accumulation.dataset <- viral_accumulation.dataset %>% tail(floor(n.infected.cells))
      gene.mtx.infected <- gene.mtx.infected[colnames(gene.mtx.infected) %in% names(viral_accumulation.dataset)]
    }
    
    message(paste0("viral RNA detected in ", n.cells.vir.RNA, " cells"))
    message(paste0("final number of infected cells: ", ncol(gene.mtx.infected)))

    # output matrices
    
    if (ncol(gene.mtx.infected) >= min.cells) {
      matrices[[paste0("infected_", d)]] <- gene.mtx.infected
      viral_accumulation[[paste0("infected_", d)]] <- viral_accumulation.dataset
    }
    if (ncol(gene.mtx.uninfected) >= min.cells) {
      matrices[[paste0("uninfected_", d)]] <- gene.mtx.uninfected
    }
  }
  
  # output integrated infected cells only
  
  matrices[["infected"]] <- matrices[grep("^infected", names(matrices))] %>% bind_cols()
  viral_accumulation.infected_cells <- viral_accumulation[grep("^infected", names(viral_accumulation))] %>% unname() %>% unlist()
  
  # reorder infected cells
  
  viral_accumulation.infected_cells <- viral_accumulation.infected_cells %>% sort()
  matrices[["infected"]] <- matrices[["infected"]][names(viral_accumulation.infected_cells)]
  
  viral_accumulation[["infected"]] <- viral_accumulation.infected_cells
  
  # remove rows of zeros
  
  matrices <- lapply(matrices, function(x) x[rowSums(x) > 0,])
  
  return(list(matrices = matrices,
              raw.matrices = raw.matrices,
              viral_accumulation = viral_accumulation))
}

#' process.gene.count.Seurat
#' 
#' Processes Seurat objects to work with the workflow, not recommended as it stores sparse matrices as full matrices
#' 
#' @param seurat.obj Seurat object
#' @param split.by.lst List to split the Seurat object
#' @param verbose Verbose
#' @param assay Seurat object assay to use
#' @param slot Seurat object slot to use
#' 
#' @return List of matrices
#' @export
#' 
#' @examples
process.gene.count.Seurat <- function(seurat.obj,
                                      split.by.lst,
                                      verbose = FALSE,
                                      assay = "RNA",
                                      slot = "counts") {
  
  seurat.lst <- list()
  matrices <- list()
  
  for (i in split.by.lst) {
    
    if (verbose) message("splitting object...")
    
    seurat.lst <- c(seurat.lst,
                    SplitObject(seurat.obj, split.by = i))
    
  }
  
  for (s in names(seurat.lst)) {
    
    if (verbose) message(paste0("processing ", s))
    
    gene.mtx <- GetAssayData(object = seurat.lst[[s]], assay = assay, slot = slot)
    
    if (ncol(gene.mtx) <= 1) next
    
    gene.mtx <- gene.mtx[rowSums(gene.mtx) > 0,] # remove rows of zeros
    
    matrices[[s]] <- gene.mtx %>% as.matrix()
  }
  return(matrices)
}

#' fit.TL
#' 
#' Fit matrices in a list to Taylor's law
#' 
#' @param matrices List of matrices
#' @param prefix.xlsx Prefix of .xlsx files for cmplxcruncher. If NULL, don't save xlsx files
#' @param zero.rate.threshold Filter out genes with percentage of zeros higher than zero.rate.threshold
#' @param normalize Boolean. Whether to transform the matrices to abundance matrix by performing cell size normalization
#' 
#' @return A list of the Taylor's parameters, the linear models and the log-transformed mean and standard deviation data
#' @export
#' 
#' @examples
fit.TL <- function(matrices,
                   prefix.xlsx = NULL,
                   zero.rate.threshold = .95,
                   normalize = TRUE) {
  
  params.tb <- data.frame(data = character(),
                          V = numeric(),
                          beta = numeric(),
                          R.squared = numeric(),
                          n.cells = numeric(),
                          n.genes = numeric(),
                          sparsity = numeric())
  log_log.mean_sd.lst <- list()
  lm.lst <- list()

  for (m in matrices %>% names()) {
    
    gene.mtx <- matrices[[m]]
    
    if (is.null(ncol(gene.mtx))) next
    
    if (!is.null(prefix.xlsx)) {
      
      # if too many cells, save three different subset files of 5000 cells
      
      if (ncol(gene.mtx) > 10000) {
        
        set.seed(12345)
        rep <- gene.mtx[sample(ncol(gene.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_1", ".xlsx"))
        
        set.seed(1000)
        rep <- gene.mtx[sample(ncol(gene.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_2", ".xlsx"))
        
        set.seed(255)
        rep <- gene.mtx[sample(ncol(gene.mtx), 5000)]
        write_xlsx(cbind(" " = rownames(rep), rep), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, "_rep_3", ".xlsx"))
        
      } else {
        write_xlsx(cbind(" " = rownames(gene.mtx), gene.mtx), paste0("cmplxcruncher_analysis/", prefix.xlsx, "_", m, ".xlsx"))
      }
    }
    if (!is.null(zero.rate.threshold)) {
      
      zero.rate.genes <- calc.0_rate(list(mtx = gene.mtx))$mtx
      zero.rate.genes <- zero.rate.genes[zero.rate.genes > zero.rate.threshold] %>% names()
      
      gene.mtx <- gene.mtx[!rownames(gene.mtx) %in% zero.rate.genes,]
    }
    if (normalize) gene.mtx <- gene.mtx %>% normalize.cells
    
    mean_sd <- data.frame(mean = rowMeans(gene.mtx),
                          sd = apply(gene.mtx, 1, sd))
    mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]

    log_log.mean_sd <- log(mean_sd)
    lm.log_log.mean_sd <- lm(sd ~ mean, data = log_log.mean_sd)
    fit.summary <- summary(lm.log_log.mean_sd)
    
    params.tb[nrow(params.tb) + 1,] <- list(m, fit.summary$coefficients[1] %>% exp(),
                                            fit.summary$coefficients[2],
                                            fit.summary$r.squared,
                                            gene.mtx %>% ncol(),
                                            log_log.mean_sd %>% nrow(),
                                            coop::sparsity(gene.mtx %>% as.matrix()))
    log_log.mean_sd.lst[[m]] <- log_log.mean_sd
    lm.lst[[m]] <- lm.log_log.mean_sd
  }
  params.tb$model <- "LLR"
  return(list(params = params.tb,
              log_log.mean_sd = log_log.mean_sd.lst,
              lm.log_log.mean_sd = lm.lst))
}

#' fit.seg.TL
#' 
#' Segmented Taylor's law fit
#' 
#' @param TL.fit.obj Object returned by fit.TL
#' @param npsi Number of breakpoints
#' @param min.cells Skip matrix if it doesn't have the minimum number of cells
#' 
#' @return A list of the Taylor's parameters, the linear models and the log-transformed mean and standard deviation data
#' @export
#' 
#' @examples
fit.seg.TL <- function(TL.fit.obj,
                       npsi = 1,
                       min.cells = 30) {
  
  params.tb <- data.frame(data = character(),
                          V = numeric(),
                          beta = numeric(),
                          R.squared = numeric(),
                          n.cells = numeric(),
                          n.genes = numeric(),
                          sparsity = numeric(),
                          model = character(),
                          breakpoint = numeric())
  log_log.mean_sd.lst <- list()
  lm.lst <- list()
  
  # get sample names
  
  datasets <- TL.fit.obj$params[[1]]
  
  for (i in datasets) {
    
    params.i <- TL.fit.obj$params[TL.fit.obj$params[1] == i & TL.fit.obj$params$model == "LLR",]

    # skip if less than min.cells cells
    
    if (params.i$n.cells < min.cells) {
      message(paste("skipped", i))
      next
    }
    
    log_log.mean_sd <- TL.fit.obj$log_log.mean_sd[[i]]
    LLR.lm.log_log.mean_sd <- lm(sd ~ mean, data = log_log.mean_sd)

    # if more than 1 breakpoint, iteratively test fit with 1 to npsi breakpoints
    
    last.lm <- LLR.lm.log_log.mean_sd
    
    for (j in 1:npsi) {
      
      curr.lm <- tryCatch(segmented(LLR.lm.log_log.mean_sd,
                                    npsi = j,
                                    control = seg.control(alpha = 0)),
                          error = function(cond) {
                            message(paste("An error ocurred while fitting segmented lm to", i, "with npsi =", j))
                            message(cond)
                            return(NULL)
                          },
                          warning = function(cond) {
                            if (cond$message == "No breakpoint estimated") {
                              message(paste("\"No breakpoint estimated\" warning while fitting segmented lm to", i, "with npsi =", j))
                              return(NULL)
                            } else {
                              message(paste("A warning was generated while fitting segmented lm to", i, "with npsi =", j))
                              message(cond)
                              #return(NULL)
                            }
                          })
      
      # check if a breakpoint was estimated
      
      if (is.null(curr.lm)) next
      
      partial.f_test0 <- anova(LLR.lm.log_log.mean_sd, curr.lm)$`Pr(>F)`[[2]]
      partial.f_test <- anova(last.lm, curr.lm)$`Pr(>F)`[[2]]
      
      # store current model if it's better than first and last model
      
      if (partial.f_test0 < .05 & partial.f_test < .05) {
        
        best.lm <- curr.lm
        best.model <- paste0("SLLR_nPSI", j)
        
        best.V <- intercept(best.lm)$mean[,1] %>% exp() %>% paste0(collapse = ",")
        best.beta <- slope(best.lm)$mean[,1] %>% paste0(collapse = ",")
      
      # if not, check whether it's the first model
        
      } else if (identical(last.lm, LLR.lm.log_log.mean_sd)) {
        
        best.lm <- LLR.lm.log_log.mean_sd
        best.model <- "LLR"
        
        best.V <- best.fit.summary$coefficients[1] %>% exp()
        best.beta <- best.fit.summary$coefficients[2]
      }
      best.fit.summary <- best.lm %>% summary()
      last.lm <- curr.lm
    }
    
    # get breakpoints
    
    breakpoints <- best.lm$psi[,2] %>% paste0(collapse = ",")
    
    params.tb[nrow(params.tb) + 1,] <- list(i, best.V,
                                            best.beta,
                                            best.fit.summary$r.squared,
                                            params.i[["n.cells"]],
                                            params.i[["n.genes"]],
                                            params.i[["sparsity"]],
                                            best.model,
                                            breakpoints)
    log_log.mean_sd.lst[[i]] <- log_log.mean_sd
    lm.lst[[i]] <- best.lm
  }
  return(list(params = params.tb,
              log_log.mean_sd = log_log.mean_sd.lst,
              lm.log_log.mean_sd = lm.lst))
}

#' pseudotime_bin
#' 
#' Bin the data into n bins, and for each one, fit to Taylor's law
#' 
#' @param gene.mtx A gene matrix
#' @param n.bins Number of bins
#' @param zero.rate.threshold Filter out genes with percentage of zeros higher than zero.rate.threshold
#' @param normalize Boolean. Whether to transform the matrices to abundance matrix by performing cell size normalization
#' 
#' @return A list of the bins, Taylor's parameters, the linear models and the log-transformed mean and standard deviation data
#' @export
#' 
#' @examples
pseudotime_bin <- function(gene.mtx,
                           n.bins = 30,
                           zero.rate.threshold = .95, # zero rate threshold
                           normalize = TRUE) {
  
  params.tb <- data.frame(pseudotime_bin = numeric(),
                          V = numeric(),
                          beta = numeric(),
                          R.squared = numeric(),
                          n.cells = numeric(),
                          n.genes = numeric(),
                          sparsity = numeric())
  
  log_log.mean_sd.lst <- list()
  lm.lst <- list()
  mtx.lst <- list()
  
  bins <- ntile(1:ncol(gene.mtx), n.bins)
  
  for (b in 1:n.bins) {
    
    gene.mtx.bin <- gene.mtx[bins == b]
    gene.mtx.bin <- gene.mtx.bin[rowSums(gene.mtx.bin) > 0,] # remove rows of zero
    
    if (!is.null(zero.rate.threshold)) {

      zero.rate.genes <- calc.0_rate(list(mtx = gene.mtx.bin))$mtx
      zero.rate.genes <- zero.rate.genes[zero.rate.genes > zero.rate.threshold] %>% names()
      
      gene.mtx.bin <- gene.mtx.bin[!rownames(gene.mtx.bin) %in% zero.rate.genes,]
      
    }
    if (normalize) gene.mtx.bin <- gene.mtx.bin %>% normalize.cells()
    
    mean_sd <- data.frame(mean = rowMeans(gene.mtx.bin),
                          sd = apply(gene.mtx.bin, 1, sd))
    mean_sd <- mean_sd[mean_sd$mean > 0 & mean_sd$sd > 0,]
    
    log_log.mean_sd <- log(mean_sd)
    lm.log_log.mean_sd <- lm(sd ~ mean, data = log_log.mean_sd)
    fit.summary <- summary(lm.log_log.mean_sd)
    
    params.tb[nrow(params.tb) + 1,] <- list(b, fit.summary$coefficients[1] %>% exp(),
                                            fit.summary$coefficients[2],
                                            fit.summary$r.squared,
                                            gene.mtx.bin %>% ncol(),
                                            log_log.mean_sd %>% nrow(),
                                            coop::sparsity(gene.mtx.bin %>% as.matrix()))
    log_log.mean_sd.lst[[b]] <- log_log.mean_sd
    lm.lst[[b]] <- lm.log_log.mean_sd
    mtx.lst[[as.character(b)]] <- gene.mtx.bin
  }
  params.tb$model <- "LLR"
  return(list(matrices = mtx.lst,
              params = params.tb,
              log_log.mean_sd = log_log.mean_sd.lst,
              lm.log_log.mean_sd = lm.lst))
}

#' calc.0_rate
#' 
#' Calculate the zero rate of genes
#' 
#' @param matrices List of gene matrices
#' 
#' @return A named list with the zero rate of each gene
#' @export
#' 
#' @examples
calc.0_rate <- function(matrices) {
  
  zero_rate.lst <- list()
  
  for (m in names(matrices)) {
    
    gene.mtx <- matrices[[m]]
    zero_rate <- rowSums(gene.mtx == 0) / ncol(gene.mtx)
    
    zero_rate.lst[[m]] <- zero_rate
  }
  return(zero_rate.lst)
}

#' plot.power_law
#' 
#' Plot Taylor's plot with metadata 
#' 
#' @param log_log.mean_sd.matrices List of means and standard deviations in log-log
#' @param genes.metadata.lst List containing, for each log-log mean-standard deviation matrix, a named list where names are the genes
#' @param legend.title Legend title
#' @param metadata.range Range on legend. If NULL, assume values on genes.metadata.lst are discrete
#' @param l10 If input is natural log, transform to log10
#' 
#' @return A ggplot
#' @export
#' 
#' @examples
plot.power_law <- function(log_log.mean_sd.matrices,
                           genes.metadata.lst,
                           legend.title = NULL,
                           metadata.range = c(0, 1),
                           l10 = TRUE) {
  
  log_log.df <- log_log.mean_sd.matrices %>% bind_rows(.id = "title")
  log_log.df$genes <- lapply(1:length(log_log.mean_sd.matrices),
                             function(x) log_log.mean_sd.matrices[[x]] %>% rownames()) %>% unlist()
  log_log.df$title <- factor(log_log.df$title, levels = names(log_log.mean_sd.matrices))
  
  value <- stack(genes.metadata.lst)
  value$genes <-  genes.metadata.lst %>% lapply(function(x) names(x)) %>% unlist()
  
  log_log.df$value <- value[match(paste0(log_log.df$title, log_log.df$genes),
                                  paste0(value$ind, value$genes)), "values"]
  
  if (l10) log_log.df <- log_log.df %>% mutate(mean = log10(exp(mean)), sd = log10(exp(sd)))
  
  gg <- ggplot(log_log.df, aes(mean, sd, color = value)) +
    facet_wrap(~title) +
    geom_point(shape = 19, stroke = 0)
  
  if (!is.null(legend.title)) gg <- gg + labs(color = legend.title)
  else gg <- gg + theme(legend.title = element_blank())
  
  if (!is.null(metadata.range)) gg <- gg + scale_color_gradient(limits = metadata.range)
  else gg <- gg + scale_color_discrete(na.translate = FALSE)
  return(gg)
}

#' process.rank
#' 
#' Convert matrices to rank matrices
#' 
#' @param cpxcruncher_corrank.path List of paths to get rank from cmplxcruncher rank file
#' @param cpxcruncher.sample_names Look for rank files containing this name
#' @param matrices List of input matrices
#' @param method Method for ties when calculating ranks
#' 
#' @return List of rank matrices
#' @export
#' 
#' @examples
process.rank <- function(cpxcruncher_corrank.path = NULL,
                         cpxcruncher.sample_names = "",
                         matrices = NULL,
                         method = "average") {
  
  rank.matrices <- list()
  
  if (!is.null(cpxcruncher_corrank.path)) {
    
    rank.files <- list.files(path = cpxcruncher_corrank.path, pattern = paste0(cpxcruncher.sample_names, ".*rank.xlsx")) %>% as.data.frame()
    
    for (f in rank.files[[1]]) {
      
      m <- f %>% gsub(".*matrix_", "", .) %>% gsub("\\.tbl.*", "", .)
      
      rank.mtx <- read_excel(paste0(cpxcruncher_corrank.path, f)) %>% as.data.frame()
      rownames(rank.mtx) <- rank.mtx[[1]]
      rank.mtx <- rank.mtx[-1]
      
      rank.matrices[[m]] <- rank.mtx
    }
  }
  
  else if(!is.null(matrices)) {
    
    for (m in names(matrices)) {
      if (method == "random") set.seed(1024)
      rank.mtx <- -matrices[[m]] %>% apply(2, rank, ties.method = method)
      rank.matrices[[m]] <- rank.mtx
    }
  }
  return(rank.matrices)
}

#' calc.RSI
#' 
#' Calculate rank stability from rank matrices
#' 
#' @param rank.matrices Input rank matrices
#' @param power_index Power index when calculating the RSI
#' @param shuffle.rep Number of replicates to calculate punctual rank stability, if NULL just calculate RSI
#' @param method Method for adjusting p-values
#' @param cdf.method Method for estimating probability density function. "ecdf" (empirical CDF) or "monoH.FC" (spline estimation)
#' @param get.rep Return replicates
#' @param keep.borders Keep first and last point when calculating punctual rank stability
#' @param genes.to.keep If not NULL, only keep these genes
#' 
#' @return A list of the RSI of each gene, or if shuffle.rep is not NULL, a table containing the RSI, PSI, mean RSI of replicates and p-values for punctual stability
#' @export
#' 
#' @examples
calc.RSI <- function(rank.matrices,
                     power_index = 4,
                     shuffle.rep = NULL,
                     method = "BH",
                     cdf.method = "ecdf",
                     get.rep = FALSE,
                     keep.borders = TRUE,
                     genes.to.keep = NULL) {
  
  RSI.lst <- list()
  
  for (m in names(rank.matrices)) {
    
    rank.mtx <- rank.matrices[[m]] %>% as.data.frame()
    
    if (!is.null(genes.to.keep)) rank.mtx <- rank.mtx[rownames(rank.mtx) %in% genes.to.keep,]
    
    max_hop <- (ncol(rank.mtx) - 1) * (nrow(rank.mtx) - 1)
    RSI <- apply(rank.mtx, 1, function(x) (1 - (sum(abs(diff(x))) / max_hop)) ** power_index)
    
    if (is.null(shuffle.rep)) {
      
      # reorder by accumulated rank (not accumulated abundance rank)
      
      stab <- apply(rank.mtx, 1, sum)
      RSI.lst[[m]] <- RSI[order(stab)]

    } else {
      
      RSI.rep <- as.data.frame(matrix(nrow = nrow(rank.mtx), ncol = 0))
      rownames(RSI.rep) <- rownames(rank.mtx)
      
      RSI.mtx <- as.data.frame(matrix(nrow = nrow(rank.mtx), ncol = 0))
      rownames(RSI.mtx) <- rownames(rank.mtx)
      
      for (i in 1:shuffle.rep) {
        
        set.seed(i)
        if (keep.borders) rank.mtx.rep <- rank.mtx[c(1, sample(2:(ncol(rank.mtx) - 1)), ncol(rank.mtx))]
        else rank.mtx.rep <- rank.mtx[sample(ncol(rank.mtx))]
        RSI.rep[i] <- apply(rank.mtx.rep, 1, function(x) (1 - (sum(abs(diff(x))) / max_hop)) ** power_index)
        
      }
      
      if (get.rep) {
        
        RSI.lst[[m]] <- RSI.rep
        
      } else {
        
        RSI.mtx$RSI <- RSI
        RSI.mtx$RSI.rep.mean <- apply(RSI.rep, 1, mean)
        RSI.mtx$RSI.boot <- RSI.mtx$RSI / RSI.mtx$RSI.rep.mean
        
        RSI.mtx$RSI.rep.sd <- apply(RSI.rep, 1, sd)
        RSI.mtx$z.score <- (RSI.mtx$RSI - RSI.mtx$RSI.rep.mean) / RSI.mtx$RSI.rep.sd
        
        RSI.mtx$p.value <- pnorm(RSI.mtx$z.score, RSI.mtx$RSI.rep.mean, RSI.mtx$RSI.rep.sd, FALSE)
        RSI.mtx$p.adjust <- p.adjust(RSI.mtx$p.value, method = method)
        
        # also add p.value based on a survival function (1 - CDF) of the RSI rep distribution , i.e., the probability of finding a RSI value at least this high
        
        RSI.mtx$p.value.survival <- lapply(1:nrow(RSI.mtx), function(x) {
          
          if (cdf.method == "ecdf") {
            
            CDF <- RSI.rep[x,] %>% unlist() %>% as.vector() %>% ecdf()
            P <- CDF(RSI.mtx$RSI[x])
            
          } else if (cdf.method == "monoH.FC") {
            
            CDF <- RSI.rep[x,] %>% unlist() %>% as.vector() %>% density() %>% with(splinefun(x, cumsum(y) / sum(y), method = "monoH.FC")) # based on https://stats.stackexchange.com/questions/78711/how-to-find-estimate-probability-density-function-from-density-function-in-r
            P <- CDF(RSI.mtx$RSI[x])
            P <- P %>% min(1) %>% max(0) # due to spline estimation, P can be < 0 or > 1
          }
          return(1 - P)
        }) %>% unlist()
        
        RSI.mtx$p.adjust.survival <- p.adjust(RSI.mtx$p.value.survival, method = method)
        
        # reorder by accumulated rank (not accumulated abundance rank)
        
        stab <- apply(rank.mtx, 1, sum)
        RSI.lst[[m]] <- RSI.mtx[order(stab),]
      }
    }
  }
  return(RSI.lst)
}

#' add.value.to.GO
#' 
#' Add median values of RSI, of Hurst analysis or whatever to the result of GO analysis. Not really meant for other data
#' 
#' @param go.tb GO results table
#' @param gene.values List of values to add to table for each cell type
#' 
#' @return A table with the median value of each cell type
#' @export
#' 
#' @examples
add.value.to.GO <- function(go.tb,
                            gene.values) {
  
  go.tb$median.hBECs <- apply(go.tb, 1, function(x) {
    genes <- strsplit(x["geneID"], "/")[[1]]
    median(gene.values$hBECs[genes])
  })
  
  go.tb$median.colon <- apply(go.tb, 1, function(x) {
    genes <- strsplit(x["geneID"], "/")[[1]]
    median(gene.values$colon[genes])
  })
  
  go.tb$median.ileum <- apply(go.tb, 1, function(x) {
    genes <- strsplit(x["geneID"], "/")[[1]]
    median(gene.values$ileum[genes])
  })
  return(go.tb)
}

#' gg.bar.GO
#' 
#' Plot GO results
#' 
#' @param go.tb GO table output
#' @param gene.values Add same values to the plot
#' @param legend Legend title
#' @param widths Widths of the heatmap and values plots
#' @param y.size axis.text.y
#' 
#' @return GO results heatmap
#' @export
#' 
#' @examples
gg.bar.GO <- function(go.tb,
                      gene.values = NULL,
                      legend = "median",
                      widths = c(2, 1),
                      y.size = 6) {
  
  if (is.null(gene.values)) {
    
    gg <- ggplot(go.tb, aes(y = `Count`, x = Description %>% factor(levels = Description %>% rev()))) +
      geom_col(aes(fill = p.adjust)) +
      coord_flip() +
      scale_fill_viridis() +
      ylab(label = "count") + xlab(label = element_blank()) +
      guides(size = "none") + guides(fill = guide_colorbar("FDR"))
    
  } else {
    
    go.tb <- add.value.to.GO(go.tb, gene.values)
    
    heat.tb <- data.frame(Description = go.tb$Description,
                          median = c(go.tb$median.hBECs,
                                     go.tb$median.colon,
                                     go.tb$median.ileum),
                          cell = c(rep("hBECs", nrow(go.tb)),
                                   rep("colon", nrow(go.tb)),
                                   rep("ileum", nrow(go.tb))))
    
    heat.tb$cell <- factor(heat.tb$cell, levels = c("hBECs", "colon", "ileum"))
    
    gg1 <- ggplot(go.tb, aes(y = `Count`, x = Description %>% factor(levels = Description %>% rev()))) +
      geom_col() +
      coord_flip() +
      scale_fill_viridis() +
      ylab(label = "count") + xlab(label = element_blank()) +
      guides(size = "none") + guides(fill = guide_colorbar("FDR")) +
      theme(axis.text.y = element_text(size = y.size))
    
    gg2 <- ggplot(heat.tb, aes(y = Description %>% factor(levels = Description %>% unique() %>% rev()),
                               x = cell, fill = median)) +
      geom_tile() + theme(axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank()) +
      scale_fill_viridis() + labs(fill = legend)

    gg <- ggarrange(gg1, gg2, nrow = 1, widths = widths)
  }
  return(gg)
}


#' down.sample
#'
#' Simulate noise in uninfected cells
#' 
#' @param mtx Uninfected cells matrix
#' @param target.counts.lst A list of the target number of cells to downsample
#' 
#' @return A downsampled matrix
#' @export
#' 
#' @examples
down.sample <- function(mtx,
                        target.counts.lst) {
  
  down.mtx <- lapply(1:ncol(mtx), function(x) {
    
    pseudotime_point <- mtx[[x]] %>% setNames(nm = rownames(mtx))

    # capture n.counts mRNA based on gene frequency
    
    n.counts <- pseudotime_point %>% sum() %>% ceiling()
    target.counts <- target.counts.lst[[x]] %>% min(n.counts)
    capture <- sample.int(n.counts, target.counts)
    genes <- findInterval(capture - 1, c(0, cumsum(pseudotime_point)), rightmost.closed = TRUE)
    obs <- names(pseudotime_point)[genes] %>% table()
    down.x <- obs[match(names(pseudotime_point), names(obs))]
    down.x[is.na(down.x)] <- 0
    return(down.x)
  }) %>% setNames(nm = colnames(mtx)) %>% bind_rows() %>% as.data.frame()
  
  rownames(down.mtx) <- rownames(mtx)
  return(down.mtx[rowSums(down.mtx) > 0,])
}

#' fit.hurst
#' 
#' Estimate H exponent for a list of matrices
#' 
#' @param matrices List of matrices
#' @param gene.lst Only estimate H for these genes
#' @param d Minimum window size for calculating H
#' 
#' @return A table with the results
fit.hurst <- function(matrices,
                      gene.lst = NULL,
                      d = 50) {

  hurst.fit <- list()
  
  if (is.null(gene.lst)) {
    gene.lst <- lapply(matrices, rownames) %>% unlist() %>% unique()
  }
  
  for (gene in gene.lst) {
    H.lst <- mapply(function(x, y) {
      
      if (!gene %in% rownames(x)) {
        message(paste0("gene ", gene, " not found in matrix ", y, ", returning NA"))
        return(list(data = y, Hrs = NA, Ht = NA, Hal = NA, He = NA))
        
      } else c(data = y, hurstexp(x[gene,] %>% unlist(), display = FALSE, d))
    }, matrices, names(matrices), SIMPLIFY = FALSE) %>% bind_rows()
    hurst.fit[[gene]] <- H.lst
  }
  return(hurst.fit)
}

#' fit.hurst_bin
#' 
#' Estimate H exponent for a list of matrices for n number of overlapping bins with a ~50% overlap
#' 
#' @param matrices List of matrices
#' @param gene.lst Only estimate H for these genes
#' @param n.bins Number of bins
#' @param d Minimum window size for calculating H
#' 
#' @return A table with the results
fit.hurst_bin <- function(matrices,
                          gene.lst = NULL,
                          n.bins = 4,
                          d = 50) {
  
  bins.lst <- lapply(matrices, function(x) ntile(1:ncol(x), n.bins + 1))

  hurst.fit <- list()
  
  if (is.null(gene.lst)) {
    gene.lst <- lapply(gene.lst, rownames) %>% unique()
  }
  
  for (gene in gene.lst) {
    H.lst <- mapply(function(x, y, z) {
      if (!gene %in% rownames(x)) {
        message(paste0("gene ", gene, " not found in matrix ", y, ", returning NA"))
        return(list(data = y, Hrs = NA, Ht = NA, Hal = NA, He = NA))
      } else {
        gene.series <- x[gene,] %>% unlist()
        return(lapply(1:(length(unique(z)) - 1),
                      function(i) c(data = y, bin = i, hurstexp(gene.series[z == i | z == i + 1], display = FALSE))))
      }
    }, matrices, names(matrices), bins.lst,
    SIMPLIFY = FALSE) %>% bind_rows()
    hurst.fit[[gene]] <- H.lst
  }
  return(hurst.fit)
}

# Langevin and Fokker-Planck model equation functions
# Not included in final manuscript but kept here since these functions work for the current data

get.mean.sd <- function(bin.TL) {
  
  bin30_mean.mtx <- bin.TL$log_log.mean_sd %>%
    lapply(function(x) tibble::rownames_to_column(x) %>% as.data.frame()) %>%
    bind_rows(.id = "pseudotime") %>%
    pivot_wider(id_cols = rowname, values_from = mean, names_from = pseudotime, values_fill = NA) %>% as.data.frame()
  
  bin30_sd.mtx <- bin.TL$log_log.mean_sd %>%
    lapply(function(x) tibble::rownames_to_column(x) %>% as.data.frame()) %>%
    bind_rows(.id = "pseudotime") %>%
    pivot_wider(id_cols = rowname, values_from = sd, names_from = pseudotime, values_fill = NA) %>% as.data.frame()
  
  bin30_sqrd.mean.mtx <- bin.TL$matrices %>%
    lapply(function(x) apply(x ** 2, 1, mean)) %>% bind_rows(.id = "pseudotime") %>% pivot_longer(-pseudotime) %>%
    pivot_wider(id_cols = name, values_from = value, names_from = pseudotime, values_fill = NA) %>% as.data.frame()
  
  rownames(bin30_mean.mtx) <- bin30_mean.mtx$rowname
  bin30_mean.mtx <- bin30_mean.mtx[-1] %>% exp()
  bin30_mean.mtx[is.na(bin30_mean.mtx)] <- 0
  
  rownames(bin30_sd.mtx) <- bin30_sd.mtx$rowname
  bin30_sd.mtx <- bin30_sd.mtx[-1] %>% exp()
  bin30_sd.mtx[is.na(bin30_sd.mtx)] <- 0
  
  rownames(bin30_sqrd.mean.mtx) <- bin30_sqrd.mean.mtx$name
  bin30_sqrd.mean.mtx <- bin30_sqrd.mean.mtx[-1]
  bin30_sqrd.mean.mtx[is.na(bin30_sqrd.mean.mtx)] <- 0
  
  return(list(mean = bin30_mean.mtx,
              sqrd.mean = bin30_sqrd.mean.mtx,
              sd = bin30_sd.mtx))
}

# mean and variance based on corresponding Fokker-Planck equation

langevin.mean <- function(alpha, beta, gamma, t) alpha * exp(1) ** (beta * exp(1) ** -(gamma * t))

langevin.sqrd.mean <- function(sigma, alpha, beta, gamma, t, tmax, V) {
  
  if (sigma == .5) {
    
    (alpha * V * exp(1) ** (2 * beta * exp(1) ** (gamma * t))) * integrate(function(x) { # limits t0 tmax
      exp(1) ** -(beta * exp(1) ** -(gamma * x))
    }, 1, tmax)$value %>% 
      return()
  }
  
  else if (sigma == 1) {
    
    (alpha * exp(1) ** ((2 * beta * exp(1) ** (gamma * t)) + ((V ** 2) * t))) %>% exp() %>%
      return()
  }
}

# variance and squared mean as a function of mean, V and t when sigma = 1 (exponential distribution)

langevin.log.variance <- function(mean, V, t) (2 * log(mean)) + log(-1 + (exp(1) ** ((V ** 2) * t)))
langevin.log.sqrd.mean <- function(mean, V, t) (2 * log(mean)) + ((V ** 2) * t)

fit.fokker_planck <- function(bin.TL,
                              gene.lst = NULL,
                              breakpoints,
                              alpha = FALSE, # whether use same alpha for mean and sd
                              V = NULL, # fix V to this value
                              estimate.V = TRUE,
                              sigma = 1, # fix sigma to this value
                              lower = c(1e-20, 1e-20, -10, -1), # lower and upper bound of nlsLM
                              upper = c(1e-1, 1e-1, 10, 1),
                              beta.lst = c(10 ** -(0:5),  -10 ** -(0:5), (10 ** -(0:5)) * 5, (-10 ** -(0:5)) * 5, 1e-100, -1e-100),
                              gamma.lst = c(10 ** -(1:5),  -10 ** -(1:5), (10 ** -(1:5)) * 5, (-10 ** -(1:5)) * 5, 1e-100, -1e-100),
                              factor = 1e-5) {
  
  mean_sd.tb <- get.mean.sd(bin.TL)
  tmax = ncol(mean_sd.tb$mean)
  
  fit.lst <- list()
  
  if (!alpha) params.tb <- data.frame(gene = character(),
                                      alpha.mean = numeric(),
                                      alpha.sqrd.mean = numeric(),
                                      beta = numeric(),
                                      gamma = numeric())
  
  else params.tb <- data.frame(gene = character(),
                               alpha = numeric(),
                               beta = numeric(),
                               gamma = numeric(),
                               V = numeric())
  
  if (is.null(gene.lst)) {
    # only get genes with all mean > 0
    gene.lst <- rownames(mean_sd.tb$mean[apply(mean_sd.tb$mean, 1, function(x) all(x != 0)),])
  }
  
  for (gene in gene.lst) {
    
    # estimate V at the gene level (each point is a pseudotime)
    
    gene.mean <- mean_sd.tb$mean[gene,] %>% unlist()
    gene.sd <- mean_sd.tb$sd[gene,] %>% unlist()
    gene.sqrd.mean <- mean_sd.tb$sqrd.mean[gene,] %>% unlist()
    
    if (estimate.V) {
      
      mean_sd <- data.frame(mean = gene.mean,
                            sd = gene.sd,
                            sqrd.mean = gene.sqrd.mean,
                            t = 1:length(gene.mean))
      
      lm.mean_sd <- tryCatch(lm(log(gene.sqrd.mean / (mean ** 2)) ~ 0 + t, data = mean_sd),
                             error = function(cond) {
                               message(paste0("An error occurred when fitting ", gene, " to log(gene.sqrd.mean / (mean ** 2)) ~ 0 + t"))
                               message(cond$message)
                               return(NULL)
                             })
      if (is.null(lm.mean_sd)) next
      
      fit.summary <- summary(lm.mean_sd)
      V <- sqrt(fit.summary$coefficients[1])
      
    }
    
    if (!sigma %in% c(.5, 1) & !is.null(sigma)) next
    
    # fit nls with a few pre-selected starting parameters
    
    if (!alpha) {
      
      # create table
      
      gene.tb <- data.frame(pseudotime_bin = c(1:tmax,
                                               1:tmax),
                            V.pseudotime = V,
                            variable = c(gene.mean,
                                         mean_sd.tb$sqrd.mean[gene,] %>% unlist()),
                            select_mean = c(rep(1, tmax),
                                            rep(0, tmax)),
                            select_sd_poisson = c(rep(0, tmax),
                                                  ifelse(sigma == .5, 1, 0) %>% rep(tmax)),
                            select_sd_exponential = c(rep(0, tmax),
                                                      ifelse(sigma == 1, 1, 0) %>% rep(tmax)))
      
      model <- variable ~ langevin.mean(alpha = alpha.mean, beta = beta, gamma = gamma, t = pseudotime_bin) * select_mean +
        langevin.sqrd.mean(alpha = alpha.sqrd.mean, beta = beta, gamma = gamma, sigma = .5, V = V.pseudotime, t = pseudotime_bin, tmax = tmax) * select_sd_poisson +
        langevin.sqrd.mean(alpha = alpha.sqrd.mean, beta = beta, gamma = gamma, sigma = 1, V = V.pseudotime, t = pseudotime_bin) * select_sd_exponential
      
      nls.models <- lapply(beta.lst, function(x) mapply(function(x, y) {
        
        tryCatch(nlsLM(model, gene.tb,
                       start = list(alpha.mean = 1e-6,
                                    alpha.sqrd.mean = 1e-5,
                                    beta = x,
                                    gamma = y),
                       lower = lower,
                       upper = upper,
                       control = nls.lm.control(maxiter = 1000, factor = factor, maxfev = 1000)),
                 error = function(cond) {
                   message(paste0("An error occurred when fitting ", gene, " to Fokker-Planck with start parameters: alpha.mean = 1e-6, alpha.sd = 1e-5, beta = ", x, ", gamma = ", y))
                   message(cond$message)
                   return(NULL)
                 })
      }, x, gamma.lst)) %>% unlist(recursive = FALSE)
      
      # compare models and get the one with the least sum of squared residuals
      
      best.model <- NULL
      
      for (curr.model in nls.models) {
        
        if (is.null(curr.model)) next
        
        if (!is.null(curr.model) & is.null(best.model)) {
          best.model <- curr.model
          next
        }
        
        if (sum(summary(curr.model)$residuals ** 2) < sum(summary(best.model)$residuals ** 2)) best.model <- curr.model
      }
      
      if (is.null(best.model)) {
        message(paste0("Warning: no model fitted for ", gene))
        next
      }
      
      nls.summary <- summary(best.model)
      fit.lst[[gene]] <- best.model
      params.tb[nrow(params.tb) + 1,] <- list(gene, nls.summary$coefficients[1,1],
                                              nls.summary$coefficients[2,1],
                                              nls.summary$coefficients[3,1],
                                              nls.summary$coefficients[4,1],
                                              V)
    } else {
      
      # create table
      # only fit the mean and later, calculate variance as a function of mean, V and t
      # here, sigma = 1
      
      gene.tb <- data.frame(pseudotime_bin = 1:tmax,
                            V.pseudotime = V,
                            variable = gene.mean)
      
      model <- variable ~ langevin.mean(alpha = alpha, beta = beta, gamma = gamma, t = pseudotime_bin)
      
      nls.models <- lapply(beta.lst, function(x) mapply(function(x, y) {
        
        tryCatch(nlsLM(model, gene.tb,
                       start = list(alpha = 1e-2,
                                    beta = x,
                                    gamma = y),
                       lower = lower,
                       upper = upper,
                       control = nls.lm.control(maxiter = 1000, factor = factor, maxfev = 1000)),
                 error = function(cond) {
                   message(paste0("An error occurred when fitting ", gene, " to Fokker-Planck with start parameters: alpha = 1e-2, beta = ", x, ", gamma = ", y))
                   message(cond$message)
                   return(NULL)
                 })
      }, x, gamma.lst)) %>% unlist(recursive = FALSE)
      
      # compare models and get the one with the least sum of squared residuals
      
      best.model <- NULL
      
      for (curr.model in nls.models) {
        
        if (is.null(curr.model)) next
        
        if (!is.null(curr.model) & is.null(best.model)) {
          best.model <- curr.model
          next
        }
        
        if (sum(summary(curr.model)$residuals ** 2) < sum(summary(best.model)$residuals ** 2)) best.model <- curr.model
      }
      
      if (is.null(best.model)) {
        message(paste0("Warning: no model fitted for ", gene))
        next
      }
      
      nls.summary <- summary(best.model)
      fit.lst[[gene]] <- best.model
      params.tb[nrow(params.tb) + 1,] <- list(gene, nls.summary$coefficients[1,1],
                                              nls.summary$coefficients[2,1],
                                              nls.summary$coefficients[3,1],
                                              V)
    }
  }
  return(list(fit = fit.lst,
              params = params.tb))
}

plot.langevin <- function(bin.TL,
                          langevin = c("pred mean", "pred SD"),
                          gene,
                          breakpoints,
                          alpha = NULL,
                          alpha.sqrd.mean = 1,
                          alpha.mean = 1,
                          beta = 1,
                          gamma = 1,
                          V = 0) {
  
  mean_sd.tb <- get.mean.sd(bin.TL)
  
  gene.mean <- mean_sd.tb$mean[gene,] %>% unlist()
  gene.sd <- mean_sd.tb$sd[gene,] %>% unlist()
  
  if (!is.null(alpha)) {
    alpha.mean <- alpha
    alpha.sqrd.mean <- alpha ** 2
  }
  
  tmax <- ncol(mean_sd.tb$mean)
  sigma <- 1
  
  gg.table <- data.frame(`obs mean` = gene.mean,
                         `obs SD` = gene.sd,
                         `pred mean` = lapply(1:tmax, function(x) {
                           
                           langevin.mean(alpha.mean, beta, gamma, x)
                           
                         }) %>% unlist(),
                         check.names = FALSE)
  
  # if using the same alpha, calculate SD as a function of mean, V and t, where sigma is assumed to be 1
  
  if (!is.null(alpha)) {
    
    log.SD <- lapply(1:tmax, function(x) {
      langevin.log.variance(gg.table$`pred mean`[[x]], V, x) %>% exp()
    }) %>% unlist()
    
    gg.table$`pred SD` <- sqrt(log.SD)
    
  } else {
    
    gg.table$`pred sqrd mean` <- lapply(1:tmax, function(x) {
      
      langevin.sqrd.mean(sigma, alpha.sqrd.mean, beta, gamma, x, tmax, V)
      
    }) %>% unlist()
    gg.table$`pred SD` <- sqrt(gg.table$`pred sqrd mean` - (gg.table$`pred mean` ** 2))
  }
  
  if (is.null(langevin)) gg.table <- gg.table[c("obs mean", "obs SD")]
  else gg.table <- gg.table[c("obs mean", "obs SD", langevin)]
  
  gg.table <- gg.table %>% tibble::rownames_to_column() %>% pivot_longer(-rowname)
  gg.table$rowname <- as.numeric(gg.table$rowname)
  
  gg <- ggplot(gg.table, aes(x = rowname, y = value, color = name)) + geom_line() +
    theme(legend.title = element_blank(), axis.title = element_blank())
  if ("pred mean" %in% langevin | "pred SD" %in% langevin) gg <- gg + 
    labs(subtitle = paste0(" alpha (mean) = ", alpha.mean,
                           " alpha (squared mean) = ", alpha.sqrd.mean,
                           " beta = ", beta,
                           " gamma = ", gamma)) +
    theme(plot.subtitle = element_text(size = 7))
  
  return(gg)
}
