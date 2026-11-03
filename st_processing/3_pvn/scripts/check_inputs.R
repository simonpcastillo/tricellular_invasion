# Load the repository configuration; set PVN_CONFIG to use a different file.
.bootstrap_path <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (!length(.bootstrap_path)) stop('Run with Rscript, not source(): cannot infer script path')
.bootstrap_path <- normalizePath(sub('^--file=', '', .bootstrap_path[1]), mustWork = TRUE)
.project_root <- normalizePath(file.path(dirname(.bootstrap_path), '..'), mustWork = TRUE)
.config_file <- Sys.getenv('PVN_CONFIG', file.path(.project_root, 'config', 'config.R'))
source(.config_file, local = FALSE)

required_packages <- c('data.table','ggplot2','patchwork','RANN','spatstat.geom',
 'scales','scattermore','svglite','ggh4x','cluster','dplyr','tidyr','stringr',
 'Seurat','SeuratObject','Matrix','readr')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly=TRUE)]
checks <- c(setNames(FILES,paste0('qupath:',names(FILES))),
            setNames(EXPR,paste0('expression:',names(EXPR))),
            region_fractions=REGION_FRACTIONS_CSV,
            setNames(file.path(HELPERS_DIR,c('pvn_enrichment_tests.R',
               'pvn_vs_global_tests.R','pvn_composition_heatmap.R')),
             c('helper:enrichment','helper:global','helper:heatmap')))
cat('PROJECT_ROOT:',PROJECT_ROOT,'\nDATA_ROOT:',DATA_ROOT,
    '\nOUT_ROOT:',ROOT_OUT,'\n\n')
for (i in seq_along(checks)) cat(if(file.exists(checks[i])) 'OK ' else 'MISSING ',
                                  names(checks)[i],': ',checks[i],'\n',sep='')
if (length(missing_packages)) cat('Missing R packages: ',
                                    paste(missing_packages,collapse=', '),'\n')
if (length(missing_packages)||any(!file.exists(checks)))
  stop('Preflight failed. Fix missing dependencies before running all steps.',call.=FALSE)
message('All mandatory inputs and packages found.')
