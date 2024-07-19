# System and transcript dynamics of cells infected with severe acute respiratory syndrome virus 2 (SARS-CoV-2)

### Description
R scripts to perform the analysis of system and transcript dynamics along the course of SARS-CoV-2 infection using scRNA-seq data. In brief, human bronchial epithelial cells (hBECs) and cells from human colon and ileum organoids infected with SARS-CoV-2 are ordered in infection pseudotime according to their viral RNA accumulation. The system dynamics along the infection is investigated by fitting the data to Taylor's power law, and the dynamics of the cellular transcripts in response to infection is investigated by their rank stability and auto-correlation.

### Overview
To run the analysis, pre-processed data (raw counts and transformed data) must be obtained from Gutiérrez and Elena 2022 (doi: 10.1038/s42003-022-04253-4), which is available at the [Zenodo repository](https://zenodo.org/records/7198900). Please, place the raw count matrices (from Suppl_Data_1) and transformed (expected) TPT matrices (from Suppl_Data_2) at the root. Please, also place the HuRI interactome (HuRI.tsv), available [here](http://www.interactome-atlas.org/download), at the root folder; and place at the prev_res folder (created automatically by source_me.R) both the Seurat object (COVID19_July.rda) from Triana et al. (doi: 10.15252/msb.202110232), available [here](https://doi.org/10.6084/m9.figshare.13703752.v1), and the supplementary file (journal.pbio.3001143.s006.xlsx) from Ravindra et al. (doi: 10.1371/journal.pbio.3001143). Below, a short explanation of the scripts in the order they should be run.

#### Contents
##### source_me.R
Some control variables and main functions. The number of bins for the analysis (default: 30) is set for each cell type in this file. As an example, some figures for the analysis with 10 and 50 bins are available at the example_figs folder.

##### fit_taylors_law.R
Segmented and unsegmented fit using all cells.

##### fit_taylors_law_pseudotime.R
Segmented and unsegmented fit using binned data.

##### gene_behaviour.R
Gene rank stability across all cells ordered according to viral RNA accumulation.

##### gene_behaviour_pseudotime.R
Gene rank stability using the means at each bin and rank auto-correlation analyses.

The analysis was performed on the following R environment:

```
> sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin22.4.0 (64-bit)
Running under: macOS Ventura 13.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /opt/homebrew/Cellar/r/4.3.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] igraph_1.5.1          gridExtra_2.3         pracma_2.4.2          ggrepel_0.9.3         viridis_0.6.4         viridisLite_0.4.2     GO.db_3.17.0         
 [8] org.Hs.eg.db_3.17.0   AnnotationDbi_1.62.2  IRanges_2.34.1        S4Vectors_0.38.2      Biobase_2.60.0        BiocGenerics_0.46.0   clusterProfiler_4.8.3
[15] SeuratObject_4.1.3    Seurat_4.3.0.1        ggpubr_0.6.0          segmented_1.6-4       nlme_3.1-163          MASS_7.3-60           readxl_1.4.3         
[22] writexl_1.4.2         ggplot2_3.4.3         tidyr_1.3.0           dplyr_1.1.3          

loaded via a namespace (and not attached):
  [1] fs_1.6.3                    matrixStats_1.0.0           spatstat.sparse_3.0-2       bitops_1.0-7                enrichplot_1.20.3          
  [6] HDO.db_0.99.1               httr_1.4.7                  RColorBrewer_1.1-3          insight_0.19.5              tools_4.3.1                
 [11] sctransform_0.4.0           backports_1.4.1             utf8_1.2.3                  R6_2.5.1                    lazyeval_0.2.2             
 [16] uwot_0.1.16                 mgcv_1.9-0                  withr_2.5.0                 sp_2.0-0                    prettyunits_1.2.0          
 [21] progressr_0.14.0            cli_3.6.1                   spatstat.explore_3.2-3      scatterpie_0.2.1            sandwich_3.0-2             
 [26] isoband_0.2.7               labeling_0.4.3              mvtnorm_1.2-3               spatstat.data_3.0-1         ggridges_0.5.4             
 [31] pbapply_1.7-2               yulab.utils_0.1.0           gson_0.1.0                  DOSE_3.26.1                 parallelly_1.36.0          
 [36] rstudioapi_0.15.0           RSQLite_2.3.1               generics_0.1.3              gridGraphics_0.5-1          ica_1.0-3                  
 [41] spatstat.random_3.1-6       car_3.1-2                   Matrix_1.6-1.1              fansi_1.0.4                 abind_1.4-5                
 [46] lifecycle_1.0.3             multcomp_1.4-25             carData_3.0-5               SummarizedExperiment_1.30.2 qvalue_2.32.0              
 [51] Rtsne_0.16                  blob_1.2.4                  promises_1.2.1              crayon_1.5.2                miniUI_0.1.1.1             
 [56] lattice_0.21-8              cowplot_1.1.1               KEGGREST_1.40.0             pillar_1.9.0                fgsea_1.26.0               
 [61] GenomicRanges_1.52.0        estimability_1.4.1          future.apply_1.11.0         codetools_0.2-19            fastmatch_1.1-4            
 [66] leiden_0.4.3                glue_1.6.2                  downloader_0.4              ggfun_0.1.3                 data.table_1.14.8          
 [71] vctrs_0.6.3                 png_0.1-8                   treeio_1.24.3               cellranger_1.1.0            gtable_0.3.4               
 [76] datawizard_0.9.0            cachem_1.0.8                S4Arrays_1.0.6              mime_0.12                   coop_0.6-3                 
 [81] tidygraph_1.2.3             coda_0.19-4                 survival_3.5-7              SingleCellExperiment_1.22.0 TH.data_1.1-2              
 [86] ellipsis_0.3.2              fitdistrplus_1.1-11         ROCR_1.0-11                 ggtree_3.8.2                bit64_4.0.5                
 [91] progress_1.2.2              RcppAnnoy_0.0.21            GenomeInfoDb_1.36.3         irlba_2.3.5.1               KernSmooth_2.23-22         
 [96] colorspace_2.1-0            DBI_1.1.3                   tidyselect_1.2.0            emmeans_1.10.0              bit_4.0.5                  
[101] compiler_4.3.1              DelayedArray_0.26.7         plotly_4.10.2               bayestestR_0.13.1           shadowtext_0.1.2           
[106] scales_1.2.1                lmtest_0.9-40               stringr_1.5.0               digest_0.6.33               goftest_1.2-3              
[111] spatstat.utils_3.0-3        XVector_0.40.0              htmltools_0.5.6             pkgconfig_2.0.3             MatrixGenerics_1.12.3      
[116] fastmap_1.1.1               rlang_1.1.1                 htmlwidgets_1.6.2           shiny_1.7.5                 farver_2.1.1               
[121] zoo_1.8-12                  jsonlite_1.8.7              BiocParallel_1.34.2         GOSemSim_2.26.1             RCurl_1.98-1.12            
[126] magrittr_2.0.3              GenomeInfoDbData_1.2.10     ggplotify_0.1.2             patchwork_1.1.3             parameters_0.21.2          
[131] munsell_0.5.0               Rcpp_1.0.11                 ape_5.7-1                   reticulate_1.32.0           stringi_1.7.12             
[136] ggraph_2.1.0                zlibbioc_1.46.0             plyr_1.8.8                  MAST_1.26.0                 parallel_4.3.1             
[141] listenv_0.9.0               deldir_1.0-9                Biostrings_2.68.1           graphlayouts_1.0.1          splines_4.3.1              
[146] tensor_1.5                  hms_1.1.3                   spatstat.geom_3.2-5         ggsignif_0.6.4              effectsize_0.8.6           
[151] reshape2_1.4.4              BiocManager_1.30.22         tweenr_2.0.2                httpuv_1.6.11               RANN_2.6.1                 
[156] purrr_1.0.2                 polyclip_1.10-4             future_1.33.0               scattermore_1.2             ggforce_0.4.1              
[161] broom_1.0.5                 xtable_1.8-4                tidytree_0.4.5              rstatix_0.7.2               later_1.3.1                
[166] tibble_3.2.1                aplot_0.2.1                 memoise_2.0.1               cluster_2.1.4               globals_0.16.2
```
