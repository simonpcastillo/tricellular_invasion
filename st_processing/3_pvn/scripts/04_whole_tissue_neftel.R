# Load the repository configuration; set PVN_CONFIG to use a different file.
.bootstrap_path <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (!length(.bootstrap_path)) stop('Run with Rscript, not source(): cannot infer script path')
.bootstrap_path <- normalizePath(sub('^--file=', '', .bootstrap_path[1]), mustWork = TRUE)
.project_root <- normalizePath(file.path(dirname(.bootstrap_path), '..'), mustWork = TRUE)
.config_file <- Sys.getenv('PVN_CONFIG', file.path(.project_root, 'config', 'config.R'))
source(.config_file, local = FALSE)

# Whole-tissue (not EC-anchored) Neftel composition by annotated Parent region.
# Self-contained: inputs and output paths are resolved from config/config.R.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(svglite); library(dplyr)
})
dir.create(OUT_FIG, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_CSV, recursive=TRUE, showWarnings=FALSE)
df0 <- data.table::rbindlist(lapply(names(FILES), function(img) {
  f <- FILES[[img]]
  if (!file.exists(f)) stop('Missing QuPath input: ',f)
  d <- data.table::fread(f)
  if (!all(c('Classification','Parent') %in% names(d)))
    stop('Whole-tissue plot requires Classification and Parent columns: ', f)
  d[,Image:=img][]
}), fill=TRUE)
NEFTEL <- c('MES_like','AC_like','OPC_like','NPC_like')
regions <- c('bulk','infiltrative','peritumoral')
counts <- df0[Classification %in% NEFTEL & Parent %in% regions,
              .N,by=.(Image,Parent,neftel_state=Classification)]
full_grid <- data.table::CJ(Image=names(FILES), Parent=regions,
                             neftel_state=NEFTEL, unique=TRUE)
counts <- counts[full_grid, on=.(Image,Parent,neftel_state)]
counts[is.na(N),N:=0L]
counts[, total_neftel:=sum(N),by=.(Image,Parent)]
counts[, proportion:=fifelse(total_neftel>0,N/total_neftel,NA_real_)]
profile <- counts[,.(proportion=if(all(is.na(proportion))) NA_real_ else
                       mean(proportion,na.rm=TRUE),
                    n_samples=sum(is.finite(proportion))),
                  by=.(Parent,neftel_state)]
profile[, Parent:=factor(Parent,levels=regions)]
profile[, neftel_state:=factor(neftel_state,
                              levels=c('AC_like','MES_like','NPC_like','OPC_like'))]
fwrite(counts,file.path(OUT_CSV,'whole_tissue_neftel_counts.csv'))
fwrite(profile,file.path(OUT_CSV,'whole_tissue_neftel_composition.csv'))
p <- ggplot(as.data.frame(profile[Parent %in% c('bulk','infiltrative')]),
            aes(Parent,proportion,fill=neftel_state)) +
  geom_col(width=.75) + scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  theme_classic(base_size=7) + labs(x=NULL,y='Malignant state composition',fill=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
svglite::svglite(file.path(OUT_FIG,'neftel_stack_Parent_noPVN.svg'),width=3,height=2.2)
print(p)
dev.off()
message('Whole-tissue Neftel composition saved to ',OUT_CSV)
