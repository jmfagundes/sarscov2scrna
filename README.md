# System and transcript dynamics of cells infected with severe acute respiratory syndrome virus 2 (SARS-CoV-2)

### Description
R scripts to perform the analysis of system and transcript dynamics along the course of SARS-CoV-2 infection using scRNA-seq data. In brief, human bronchial epithelial cells (hBECs) and cells from human colon and ileum organoids infected with SARS-CoV-2 are ordered in infection pseudotime according to their viral RNA accumulation. The system dynamics along the infection is investigated by fitting the data to Taylor's power law, and the dynamics of the cellular transcripts in response to infection is investigated by their rank stability and auto-correlation.

### Overview
To run the analysis, pre-processed data (raw counts and transformed data) must be obtained from Gutiérrez and Elena 2022 (doi: 10.1038/s42003-022-04253-4), which is available at the [Zenodo repository](https://zenodo.org/records/7198900). Please, place the raw count matrices (from Suppl_Data_1) and transformed (expected) TPT matrices (from Suppl_Data_2) at the root. Below, a short explanation of the scripts in the order they should be run.

#### Contents
##### source_me.R
Some control variables and main functions.

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
 [1] gridExtra_2.3         pracma_2.4.2          ggrepel_0.9.3         viridis_0.6.4         viridisLite_0.4.2     GO.db_3.17.0          org.Hs.eg.db_3.17.0  
 [8] AnnotationDbi_1.62.2  IRanges_2.34.1        S4Vectors_0.38.2      Biobase_2.60.0        BiocGenerics_0.46.0   clusterProfiler_4.8.3 SeuratObject_4.1.3   
[15] Seurat_4.3.0.1        ggpubr_0.6.0          segmented_1.6-4       nlme_3.1-163          MASS_7.3-60           readxl_1.4.3          writexl_1.4.2        
[22] ggplot2_3.4.3         tidyr_1.3.0           dplyr_1.1.3          

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.21        splines_4.3.1           later_1.3.1             ggplotify_0.1.2         bitops_1.0-7            tibble_3.2.1           
  [7] cellranger_1.1.0        datawizard_0.9.0        polyclip_1.10-4         lifecycle_1.0.3         rstatix_0.7.2           globals_0.16.2         
 [13] lattice_0.21-8          insight_0.19.5          backports_1.4.1         magrittr_2.0.3          rmarkdown_2.25          plotly_4.10.2          
 [19] httpuv_1.6.11           sctransform_0.4.0       sp_2.0-0                spatstat.sparse_3.0-2   reticulate_1.32.0       cowplot_1.1.1          
 [25] pbapply_1.7-2           DBI_1.1.3               RColorBrewer_1.1-3      multcomp_1.4-25         abind_1.4-5             zlibbioc_1.46.0        
 [31] Rtsne_0.16              purrr_1.0.2             ggraph_2.1.0            RCurl_1.98-1.12         TH.data_1.1-2           yulab.utils_0.1.0      
 [37] sandwich_3.0-2          tweenr_2.0.2            GenomeInfoDbData_1.2.10 enrichplot_1.20.3       irlba_2.3.5.1           listenv_0.9.0          
 [43] spatstat.utils_3.0-3    tidytree_0.4.5          goftest_1.2-3           spatstat.random_3.1-6   fitdistrplus_1.1-11     parallelly_1.36.0      
 [49] leiden_0.4.3            codetools_0.2-19        ggforce_0.4.1           DOSE_3.26.1             tidyselect_1.2.0        aplot_0.2.1            
 [55] farver_2.1.1            effectsize_0.8.6        matrixStats_1.0.0       spatstat.explore_3.2-3  jsonlite_1.8.7          tidygraph_1.2.3        
 [61] ellipsis_0.3.2          progressr_0.14.0        emmeans_1.10.0          ggridges_0.5.4          survival_3.5-7          tools_4.3.1            
 [67] treeio_1.24.3           ica_1.0-3               Rcpp_1.0.11             glue_1.6.2              xfun_0.40               qvalue_2.32.0          
 [73] GenomeInfoDb_1.36.3     withr_2.5.0             fastmap_1.1.1           fansi_1.0.4             digest_0.6.33           estimability_1.4.1     
 [79] gridGraphics_0.5-1      R6_2.5.1                mime_0.12               colorspace_2.1-0        scattermore_1.2         tensor_1.5             
 [85] spatstat.data_3.0-1     RSQLite_2.3.1           utf8_1.2.3              generics_0.1.3          data.table_1.14.8       graphlayouts_1.0.1     
 [91] httr_1.4.7              htmlwidgets_1.6.2       parameters_0.21.2       scatterpie_0.2.1        uwot_0.1.16             pkgconfig_2.0.3        
 [97] gtable_0.3.4            blob_1.2.4              lmtest_0.9-40           XVector_0.40.0          shadowtext_0.1.2        htmltools_0.5.6        
[103] carData_3.0-5           fgsea_1.26.0            scales_1.2.1            png_0.1-8               ggfun_0.1.3             knitr_1.44             
[109] rstudioapi_0.15.0       reshape2_1.4.4          coda_0.19-4             zoo_1.8-12              cachem_1.0.8            stringr_1.5.0          
[115] KernSmooth_2.23-22      parallel_4.3.1          miniUI_0.1.1.1          HDO.db_0.99.1           pillar_1.9.0            vctrs_0.6.3            
[121] RANN_2.6.1              promises_1.2.1          car_3.1-2               xtable_1.8-4            cluster_2.1.4           evaluate_0.21          
[127] isoband_0.2.7           mvtnorm_1.2-3           cli_3.6.1               compiler_4.3.1          rlang_1.1.1             crayon_1.5.2           
[133] future.apply_1.11.0     ggsignif_0.6.4          labeling_0.4.3          fs_1.6.3                plyr_1.8.8              stringi_1.7.12         
[139] deldir_1.0-9            BiocParallel_1.34.2     munsell_0.5.0           Biostrings_2.68.1       coop_0.6-3              lazyeval_0.2.2         
[145] spatstat.geom_3.2-5     bayestestR_0.13.1       GOSemSim_2.26.1         Matrix_1.6-1.1          patchwork_1.1.3         bit64_4.0.5            
[151] future_1.33.0           KEGGREST_1.40.0         shiny_1.7.5             ROCR_1.0-11             igraph_1.5.1            broom_1.0.5            
[157] memoise_2.0.1           ggtree_3.8.2            fastmatch_1.1-4         bit_4.0.5               downloader_0.4          gson_0.1.0             
[163] ape_5.7-1
```
