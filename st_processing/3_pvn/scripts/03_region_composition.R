# Load the repository configuration; set PVN_CONFIG to use a different file.
.bootstrap_path <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (!length(.bootstrap_path)) stop('Run with Rscript, not source(): cannot infer script path')
.bootstrap_path <- normalizePath(sub('^--file=', '', .bootstrap_path[1]), mustWork = TRUE)
.project_root <- normalizePath(file.path(dirname(.bootstrap_path), '..'), mustWork = TRUE)
.config_file <- Sys.getenv('PVN_CONFIG', file.path(.project_root, 'config', 'config.R'))
source(.config_file, local = FALSE)

# Independent region composition summary (formerly panel 2C in 2_hd_to_pvn.R).
# Figures 2B and 2D originated in QuPath and are not reproduced here.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(svglite)
})
if (!file.exists(REGION_FRACTIONS_CSV)) stop('Missing region fractions: ', REGION_FRACTIONS_CSV)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_CSV, recursive = TRUE, showWarnings = FALSE)
wide <- readr::read_csv(REGION_FRACTIONS_CSV, show_col_types = FALSE)
cell_cols <- c('tumour','endothelial','mg','bam','macper','tcell','neurons')
required <- c('region','sample',cell_cols)
missing <- setdiff(required, names(wide))
if (length(missing)) stop('Region fractions CSV lacks: ', paste(missing, collapse=', '))
region_levels <- c('core','im','pt')
region_labels <- c('Core','IM','PT')
region_full <- c('Core','Infiltrative margin','Peritumoral')
long <- wide |>
  tidyr::pivot_longer(dplyr::all_of(cell_cols), names_to='cell_type', values_to='fraction') |>
  dplyr::mutate(region=factor(region, levels=region_levels, labels=region_labels),
                sample=factor(sample))
if (anyNA(long$region)) warning('Unknown anatomical region labels found')
facet_order <- long |>
  dplyr::filter(region == 'Core') |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(core=mean(fraction, na.rm=TRUE), .groups='drop') |>
  dplyr::arrange(dplyr::desc(core)) |>
  dplyr::pull(cell_type)
long <- dplyr::mutate(long, cell_type=factor(cell_type, levels=facet_order))
summ <- long |>
  dplyr::group_by(cell_type,region) |>
  dplyr::summarise(n=dplyr::n(), mean=mean(fraction,na.rm=TRUE),
                    sd=sd(fraction,na.rm=TRUE), .groups='drop')
pal <- c(Core='#000000', IM='#777777',PT='#c9c9c9')
facet_labels <- c(tumour='Tumour', endothelial='Endothelial', mg='Microglia',
                  bam='BAM', macper='Periph. mac.', tcell='T cell', neurons='Neurons')
p <- ggplot(summ, aes(region,mean,fill=region)) +
  geom_col(width=0.68, colour='grey25', linewidth=0.25) +
  geom_point(data=long, aes(y=fraction,shape=sample),
             position=position_jitter(width=0.11,height=0,seed=1),
             size=1.6, stroke=0.45, fill='white',colour='grey15') +
  facet_wrap(~cell_type,nrow=1,labeller=labeller(cell_type=facet_labels)) +
  scale_fill_manual(values=pal,labels=region_full,name='Region') +
  scale_shape_manual(values=c(21,24),name='Mouse') +
  scale_y_continuous(limits=c(0,0.85),breaks=seq(0,0.8,0.2),
                     expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL,y='Fraction of cell type') +
  theme_classic(base_size=8) +
  theme(strip.background=element_blank(), strip.text=element_text(size=7,face='bold'),
        axis.text=element_text(size=6),legend.title=element_text(size=7),
        legend.text=element_text(size=7),panel.spacing=unit(0.5,'lines'))
svglite::svglite(file.path(OUT_FIG,'2C_region_fractions_faceted.svg'),width=7,height=2.5)
print(p)
dev.off()
ggsave(file.path(OUT_FIG,'2C_region_fractions_faceted.png'),p,
       width=7,height=2.5,dpi=300,bg='white')
ggsave(file.path(OUT_FIG,'2C_region_fractions_faceted.pdf'),p,
       width=7,height=2.5)
readr::write_csv(summ,file.path(OUT_CSV,'region_fractions_summary.csv'))
message('Region summary written to ', OUT_CSV)
