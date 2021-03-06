#!/usr/bin/env Rscript

library(methods)

options(ENTREZ_KEY = "a65e721e6610fc345921a5a2778178171707")

needed.packages = c('argparser', 'ggplot2', 'plotly', 'magrittr', 'dplyr', 'rstudioapi', 'taxize', 'phangorn', 'htmlwidgets', 'BiocManager', 'stringr')
absent.packages = needed.packages[!(needed.packages %in% rownames(installed.packages()))]
install.packages(absent.packages, repos = "http://cran.us.r-project.org")

needed.packages.bio = c('ggtree', "ape", "Biostrings")
absent.packages.bio = needed.packages.bio[!(needed.packages.bio %in% rownames(installed.packages()))]
if(length(absent.packages.bio) > 0) BiocManager::install(absent.packages.bio)

for(pack in needed.packages) eval(parse(text = sprintf('library(%s)', pack)))
for(pack in needed.packages.bio) eval(parse(text = sprintf('library(%s)', pack)))


### defining startup argements
p <- arg_parser("RTrans transcriptome analysis pipeline")
p <- add_argument(p, "--input-csv", short = '-i', help="Input RunInfo CSV file from NCBI SRA", default = NA)
p <- add_argument(p, "--clade-name", short = '-n', help="Name of clade (for naming plots). For example, insecta", default = 'Default')
p <- add_argument(p, "--out-dir", short = '-o', help="Output directory", default = NA)
p = parse_args(p, argv = commandArgs(trailingOnly = TRUE))


if(is.na(p$input_csv))  stop('Input csv file should be specified')
if(!file.exists(p$input_csv))  stop(sprintf('Input csv file %s is absent', p$input_csv))


if(is.na(p$out_dir)){
  p$out_dir = tryCatch(expr = {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    }, error = function (err){
      getwd()
    }, warning = function (err){
      dirname(rstudioapi::getActiveDocumentContext()$path)
    })
}
dir.create(p$out_dir, showWarnings = FALSE)

# p$input_csv = 'insecta.csv'
# p$clade_name = 'Insecta'
run.info.file.name = p$input_csv
descr = p$clade_name

run.info.table = read.delim2(run.info.file.name, sep = ',', quote = '\"', check.names = FALSE, as.is = TRUE) %>% as.data.frame
setwd(p$out_dir)

run.info.table$ReleaseDate = as.Date(run.info.table$ReleaseDate)


found.taxons = levels(as.factor(run.info.table$TaxID))
cat(sprintf('Total %d runs found in file "%s". %d NCBI taxons present\n', nrow(run.info.table), run.info.file.name, length(found.taxons)))
cat(sprintf('   %d BioSamples, %d BioProjects',
            length(levels(as.factor(run.info.table$BioSample))),
            length(levels(as.factor(run.info.table$BioProject)))
            ))

cat(sprintf('   %.2f Gbases, %.2f billion reads\n',
            sum(run.info.table$bases) / 1e+9,
            sum(run.info.table$spots) / 1e+9
))

for(c_LibrarySource in names(table(run.info.table$LibrarySource))){
  keep = run.info.table$LibrarySource %in% c_LibrarySource
  cat(sprintf('%s: Total %d runs; %.2f Gbases (%.1f%%), %.2f billion reads \n',
              c_LibrarySource,
              sum(keep),
              sum(run.info.table$bases[keep]) / 1e+9,
              sum(run.info.table$bases[keep]) / sum(run.info.table$bases) * 100,
              sum(run.info.table$spots[keep]) / 1e+9))
}


###########################################################################################
###########################################################################################

max_log10 = ceiling(log10(sum(run.info.table$bases)))
c_breaks = 10^(3:(max_log10 - 3))

gx = ggplot(run.info.table, aes(x = ReleaseDate, y = bases)) + 
  geom_hline(yintercept = c_breaks, linetype="dashed", color = "grey", size=0.5) +
  geom_point(alpha = 0.05, fill = '#000000') +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "����� ����������� ����.", trans = 'log10', breaks = c_breaks) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")

file.name = sprintf('%s - timeline simple.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, gx, device="png", width = 12, height = 8, scale = 1.0)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)

###########################################################################################
###########################################################################################

# The direction argument allows to reverse the palette
gx.kd = ggplot(run.info.table, aes(x = ReleaseDate, y = bases)) + 
  stat_density_2d(aes(fill = ..density.. ** 0.33), geom = "raster", contour = FALSE) +
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ^ 2, fill = Label)) + 
  scale_fill_distiller(palette='Spectral', direction=-1) +
  #scale_fill_gradient2(low = "white", mid = "#9ae002", high = "#e05702", midpoint = 0.04) +
  geom_hline(yintercept = c_breaks, linetype="dashed", color = "grey", size=0.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position='none') +
  scale_y_continuous(name = "����� ����������� ����.", trans = 'log10', breaks = c_breaks) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")


file.name = sprintf('%s - timeline density.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, gx.kd, device="png", width = 12, height = 8, scale = 0.7)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)

###########################################################################################
###########################################################################################


#'ReleaseDate'
df = cbind(run.info.table[,c('bases'), drop = FALSE], as.data.frame(rep(1, nrow(run.info.table))))
colnames(df) = c('bases', 'runs_count')
df.aggr = aggregate.data.frame(x = df, by = list(BioProject = run.info.table$BioProject),  FUN = function(x){ sum (x) })

all.ReleaseDates = run.info.table$ReleaseDate
names(all.ReleaseDates) = run.info.table$BioProject
all.ReleaseDates.nr = all.ReleaseDates[!duplicated(names(all.ReleaseDates))]

# all.LibTypes = run.info.table$ReleaseDate
# names(all.ReleaseDates) = run.info.table$BioProject
# all.ReleaseDates.nr = all.ReleaseDates[!duplicated(names(all.ReleaseDates))]

all.ReleaseDates.nr.df = as.data.frame(all.ReleaseDates.nr[df.aggr$BioProject])
rownames(all.ReleaseDates.nr.df) = df.aggr$BioProject
rownames(df.aggr) = df.aggr$BioProject
df.aggr.ext = cbind(df.aggr, all.ReleaseDates.nr.df)
colnames(df.aggr.ext) = c('BioProject', 'bases', 'runs_count', 'ReleaseDate')
df.aggr.ext$ReleaseDate = as.Date(df.aggr.ext$ReleaseDate)
df.aggr.ext$log10_runs_count = log10(df.aggr.ext$runs_count)

df.aggr.ext$BioProject[which.max(df.aggr.ext$bases)]

max_log10 = ceiling(log10(sum(run.info.table$bases)))
c_breaks = 10^(5:(max_log10 - 1))
gx = ggplot(df.aggr.ext, aes(x = ReleaseDate, y = bases)) + 
  geom_hline(yintercept = c_breaks, linetype="dashed", color = "grey", size=0.5) +
  geom_point(alpha = 0.15, fill = '#000000', aes(size = log10_runs_count)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "����� ����������� ����.", trans = 'log10', breaks = c_breaks, limits = c(1e+6, 2e+14)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")

file.name = sprintf('%s - timeline2 - BioProjects.png', descr)
file.name.HTML = sprintf('%s - timeline2 - BioProjects.html', descr)
ggsave(file.name, gx, device="png", width = 12, height = 8, scale = 0.7)
saveWidget(as_widget(ggplotly(gx)), file.name.HTML, selfcontained = TRUE)


#########################################################################################
##########################################################################################


#'ReleaseDate'
df = cbind(run.info.table[,c('bases'), drop = FALSE], as.data.frame(rep(1, nrow(run.info.table))))
colnames(df) = c('bases', 'runs_count')
df.aggr = aggregate.data.frame(x = df, by = list(BioProject.and.LibrarySource = sprintf('%s---%s', run.info.table$BioProject, run.info.table$LibrarySource)),  FUN = function(x){ sum (x) })

df.aggr.part = t(as.data.frame(strsplit(df.aggr$BioProject.and.LibrarySource, '---')))
colnames(df.aggr.part) = c('BioProject', 'LibrarySource')

df.aggr = cbind(df.aggr, df.aggr.part)
rownames(df.aggr) = df.aggr$BioProject.and.LibrarySource

all.ReleaseDates = run.info.table$ReleaseDate
names(all.ReleaseDates) = run.info.table$BioProject
all.ReleaseDates.nr = all.ReleaseDates[!duplicated(names(all.ReleaseDates))]

all.ReleaseDates.nr.df = as.data.frame(all.ReleaseDates.nr[df.aggr$BioProject])
rownames(all.ReleaseDates.nr.df) = df.aggr$BioProject.and.LibrarySource
df.aggr.ext = cbind(df.aggr, all.ReleaseDates.nr.df)
colnames(df.aggr.ext) = c('BioProject.and.LibrarySource', 'bases', 'runs_count', 'BioProject', 'LibrarySource', 'ReleaseDate')

df.aggr.ext$ReleaseDate = as.Date(df.aggr.ext$ReleaseDate)
df.aggr.ext$log10_runs_count = log10(df.aggr.ext$runs_count)

df.aggr.ext$BioProject[which.max(df.aggr.ext$bases)]

max_log10 = ceiling(log10(sum(run.info.table$bases)))
c_breaks = 10^(5:(max_log10 - 1))
gx = ggplot(df.aggr.ext, aes(x = ReleaseDate, y = bases)) + 
  geom_hline(yintercept = c_breaks, linetype="dashed", color = "grey", size=0.5) +
  geom_point(alpha = 0.1, aes(size = log10_runs_count, color = LibrarySource)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Total bases", trans = 'log10', breaks = c_breaks, limits = c(1e+6, 2e+14)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")

file.name = sprintf('%s - timeline2 - BioProjects ext.png', descr)
file.name.HTML = sprintf('%s - timeline2 - BioProjects ext.html', descr)
ggsave(file.name, gx, device="png", width = 12, height = 8, scale = 0.7)
saveWidget(as_widget(ggplotly(gx)), file.name.HTML, selfcontained = TRUE)


###########################################################################################
###########################################################################################





pallete = c("#F46D43", "#66C2A5", "#cd8845", "#3288BD", "#a8bf32", "#5E4FA2", "#D53E4F", "#d6d639", "#8ed384", "#9E0142", "#ebba2f")
density.cols = colorRampPalette(pallete)(30)

# 
# 
# 
# g1 = ggplot(run.info.table, aes(x = Platform, fill = LibrarySelection)) + 
#   geom_hline(yintercept = 0, linetype="dashed", color = "grey", size=0.5) +
#   geom_histogram(stat="count", aes(group = LibrarySelection)) +
#   theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
#         text = element_text(size = 12),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black")) +
#   ggtitle(sprintf('������������� ������ NCBI SRA �� ���������� - %s', descr)) + 
#   scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)")
# 
# file.name = sprintf('%s - platform distrib.png', descr)
# file.name.HTML = sprintf('%s - platform distrib.html', descr)
# ggsave(file.name, g1, device="png", width = 12, height = 8)
# saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)
# 

g2 = ggplot(run.info.table, aes(x = Platform)) + 
  geom_hline(yintercept = c(0, 1, 10, 100, 1000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  scale_fill_brewer(palette = 'Set1') +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���������� (log-�������) - %s', descr)) + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)", trans = 'log10')

file.name = sprintf('%s - platform distrib simple.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g2, device="png", width = 12, height = 8, scale = 0.7)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)

#run.info.table$LibraryStrategy

g3 = ggplot(run.info.table, aes(x = LibrarySelection)) + 
  geom_hline(yintercept = c(0, 1, 10, 100, 1000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���� ��������� (log-�������) - %s', descr)) + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)", trans = 'identity')

file.name = sprintf('%s - LibrarySelection distrib simple.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g3, device="png", width = 12, height = 8, scale = 0.7)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)

run.info.table$LibraryStrategy[which(run.info.table$LibraryStrategy %in% "Tethered Chromatin Conformation Capture")] = 'HiC'

g5 = ggplot(run.info.table, aes(x = LibraryStrategy, fill = Platform)) + 
  #geom_hline(yintercept = seq(0, 120000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���� ��������� - %s', descr)) + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)", trans = 'identity')

file.name = sprintf('%s - LibraryStrategy distrib.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g5, device="png", width = 12, height = 8, scale = 1.0)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)

run.info.table.1 = run.info.table[run.info.table$LibrarySource %in% 'METAGENOMIC',]

g5.1 = ggplot(run.info.table.1, aes(x = LibraryStrategy, fill = Platform)) + 
  #geom_hline(yintercept = seq(0, 120000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���� ��������� - %s', descr), subtitle = '������ ������������ ��������') + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*).", trans = 'identity')

file.name = sprintf('%s - LibraryStrategy distrib metagenomic only.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g5.1, device="png", width = 8, height = 6, scale = 1.0)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)






Illumina.count = sum(run.info.table$Platform %in% 'ILLUMINA')

ONT.count = sum(run.info.table$Platform %in% 'OXFORD_NANOPORE')
ONT.ratio = sum(run.info.table$Platform %in% 'OXFORD_NANOPORE') / nrow(run.info.table)
ONT.bases.ratio = sum(run.info.table$bases[run.info.table$Platform %in% 'OXFORD_NANOPORE']) / sum(run.info.table$bases)

PacBio.count = sum(run.info.table$Platform %in% 'PACBIO_SMRT')
PacBio.ratio = sum(run.info.table$Platform %in% 'PACBIO_SMRT') / nrow(run.info.table)
PacBio.bases.ratio = sum(run.info.table$bases[run.info.table$Platform %in% 'PACBIO_SMRT']) / sum(run.info.table$bases)

run.info.table..Illumina = run.info.table[run.info.table$Platform %in% 'ILLUMINA',]
ratio = nrow(run.info.table..Illumina) / nrow(run.info.table)
print(ratio)

lib.types = table(run.info.table..Illumina$LibraryStrategy)
keep = lib.types >= 0.01 * nrow(run.info.table..Illumina)
passed.lib.types = lib.types[keep]
sum(keep) / length(keep)

keep.table = run.info.table..Illumina$LibraryStrategy %in% names(lib.types[keep])
run.info.table..Illumina = run.info.table..Illumina[keep.table,]

colors.v1 <- c("#771C19", "#AA3929", "#E25033", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#000000")

colors.v2 <- c("#ff7f00",
                         "#ffff99", "#6a3d9a",
                         "#b15928", "#cab2d6",
                         "#B2DF8A", "#33A02C",
                         "#fdbf6f", "#A6CEE3",
                         "#1F78B4", "#fb9a99",
                         "#e31a1c")
colors.v3 <- c(
    "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
    "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

c_colors = colorRampPalette(colors.v2)( length(levels(as.factor(run.info.table..Illumina$Model))))

g6 = ggplot(run.info.table..Illumina, aes(x = LibraryStrategy, fill = Model)) + 
  #geom_hline(yintercept = seq(0, 120000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���� ��������� - %s', descr), subtitle = '(� ����� > 1% ���������, ������ Illumina)') + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)", trans = 'identity') + 
  scale_fill_manual(values = c_colors)

file.name = sprintf('%s - Illumina-Model distrib.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g6, device="png", width = 12, height = 8, scale = .8)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)


run.info.table..Illumina.1 = run.info.table..Illumina[run.info.table..Illumina$LibrarySource %in% 'METAGENOMIC',]

g6.1 = ggplot(run.info.table..Illumina.1, aes(x = LibraryStrategy, fill = Model)) + 
  #geom_hline(yintercept = seq(0, 120000, 10000), linetype="dashed", color = "grey", size=0.5) +
  geom_histogram(stat="count") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('������������� ������ NCBI SRA �� ���� ��������� - %s', descr), subtitle = '(� ����� > 1% ���������, ������ Illumina, ������������ ������)') + 
  scale_y_continuous(name = "����� �������� (SRR*, ERR*, DRR*)", trans = 'identity') + 
  scale_fill_manual(values = c_colors)

file.name = sprintf('%s - Illumina-Model distrib metagenomic only.png', descr)
#file.name.HTML = sprintf('%s - platform distrib.html', descr)
ggsave(file.name, g6.1, device="png", width = 8, height = 6, scale = .8)
#saveWidget(as_widget(ggplotly(g1)), file.name.HTML, selfcontained = TRUE)


################################################

libstrat.table = table(run.info.table$LibraryStrategy)
keep = libstrat.table >= 0.01 * sum(libstrat.table)
sum(keep) / length(keep)
passed.libstats = names(libstrat.table)[keep]

run.info.table..PL = run.info.table[run.info.table$LibraryStrategy %in% passed.libstats,, drop = FALSE]
ratio = nrow(run.info.table..PL) / nrow(run.info.table)
print(ratio)


c_colors = colorRampPalette(colors.v1)( length(passed.libstats))

g7 = ggplot(run.info.table..PL, aes(x = ReleaseDate, fill = NULL, color = LibraryStrategy)) + 
  geom_density(bw = 300, size = 1.0, alpha = 1.0) + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1.0, vjust = 1.0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggtitle(sprintf('��������� �����. ������ NCBI SRA �� ������� ������� - %s', descr), subtitle = '(� ����������� �� ���� ���������, � ����� > 1% ���������)') + 
  scale_y_continuous(name = "��������� �������������", trans = 'identity') + 
  scale_color_manual(values = colors.v3)

file.name = sprintf('%s - LibraryStrategy - timeline distrib.png', descr)
file.name.HTML = sprintf('%s - LibraryStrategy - timeline distrib.html', descr)
ggsave(file.name, g7, device="png", width = 12, height = 8, scale = 0.75)
saveWidget(as_widget(ggplotly(g7)), file.name.HTML, selfcontained = TRUE)


#######################################################################################
#######################################################################################
#######################################################################################


top.taxons.size = 50
taxon.counts = table(run.info.table$TaxID) %>% .[rev(order(.))]
top.taxons = names(taxon.counts)[1:top.taxons.size]
top.taxons.ratio = sum(taxon.counts[1:top.taxons.size]) / sum(taxon.counts)
cat(sprintf('Top %d taxons correspond for %.1f%% percents of datasets', top.taxons.size, top.taxons.ratio * 100))
#found.taxons = levels(as.factor(run.info.table$TaxID))

# res = ncbi_get_taxon_summary(top.taxons, key = 'a65e721e6610fc345921a5a2778178171707')
res2 = classification(top.taxons, db = 'ncbi')
tr <- class2tree(res2)

x = 4
for(x in 1:top.taxons.size){
  #tr$phylo$tip.label[x]
  keep = run.info.table$TaxID %in% names(res2)[x]
  l.keep.RNAseq = run.info.table$LibrarySource %in% c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC")
  l.keep.DNAseq = run.info.table$LibrarySource %in% c("GENOMIC", "METAGENOMIC")
  
  records.count = sum(keep)
  records.count.RNAseq = sum(keep & l.keep.RNAseq)
  records.count.DNAseq = sum(keep & l.keep.DNAseq)
  #Gbases.count = sum(run.info.table$bases[keep]) / 1e+9
  Tbases.count = sum(run.info.table$bases[keep]) / 1e+12
  Tbases.count.RNAseq = sum(run.info.table$bases[keep & l.keep.RNAseq]) / 1e+12
  Tbases.count.DNAseq = sum(run.info.table$bases[keep & l.keep.DNAseq]) / 1e+12
  tr$phylo$tip.label[x] = sprintf('%s  (%.d;  %.1f Gb; r = %.2f)', tr$phylo$tip.label[x],
                                  records.count, 1000*Tbases.count, Tbases.count.RNAseq / (Tbases.count.RNAseq + Tbases.count.DNAseq))

}

print(tr$phylo$tip.label)

png(filename = sprintf('%s - top %d taxons (by SRRs count) tree.png', descr, top.taxons.size), 
    width = 10, height = 15, units = 'cm', res = 600, pointsize = 5.0)
plot(tr$phylo)
dev.off()


#############################################################

top.taxons.size = 50

bases.per.taxon = c()

sdf = run.info.table[,c('TaxID', 'bases')]
bases.counts.df = aggregate.data.frame(run.info.table[,'bases',drop = FALSE], by = list(TaxID = run.info.table$TaxID), function(x) sum(x) )

taxon.bases.counts = bases.counts.df$bases
names(taxon.bases.counts) = bases.counts.df$TaxID
taxon.bases.counts %<>% .[rev(order(.))]

top.taxons = names(taxon.bases.counts)[1:top.taxons.size]
top.taxons.ratio = sum(taxon.bases.counts[1:top.taxons.size]) / sum(taxon.bases.counts)
cat(sprintf('Top %d taxons correspond for %.1f%% percents of bases', top.taxons.size, top.taxons.ratio * 100))
#found.taxons = levels(as.factor(run.info.table$TaxID))

# res = ncbi_get_taxon_summary(top.taxons, key = 'a65e721e6610fc345921a5a2778178171707')
res2 = classification(top.taxons, db = 'ncbi')
tr <- class2tree(res2)

x = 4
for(x in 1:top.taxons.size){
  #tr$phylo$tip.label[x]
  keep = run.info.table$TaxID %in% names(res2)[x]
  l.keep.RNAseq = run.info.table$LibrarySource %in% c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC")
  l.keep.DNAseq = run.info.table$LibrarySource %in% c("GENOMIC", "METAGENOMIC")
  
  records.count = sum(keep)
  records.count.RNAseq = sum(keep & l.keep.RNAseq)
  records.count.DNAseq = sum(keep & l.keep.DNAseq)
  #Gbases.count = sum(run.info.table$bases[keep]) / 1e+9
  Tbases.count = sum(run.info.table$bases[keep]) / 1e+12
  Tbases.count.RNAseq = sum(run.info.table$bases[keep & l.keep.RNAseq]) / 1e+12
  Tbases.count.DNAseq = sum(run.info.table$bases[keep & l.keep.DNAseq]) / 1e+12
  tr$phylo$tip.label[x] = sprintf('%s  (%.d;  %.1f Gb; r = %.2f)', tr$phylo$tip.label[x],
                                  records.count, 1000*Tbases.count, Tbases.count.RNAseq / (Tbases.count.RNAseq + Tbases.count.DNAseq))
}

print(tr$phylo$tip.label)

png(filename = sprintf('%s - top %d taxons (by SRRs count) tree.png', descr, top.taxons.size), 
    width = 10, height = 15, units = 'cm', res = 600, pointsize = 5.0)
plot(tr$phylo)
dev.off()



######################################################################################################
######################################################################################################
######################################################################################################

keep = run.info.table$Platform %in% c('OXFORD_NANOPORE', 'PACBIO_SMRT')
#keep = rep(TRUE, nrow(run.info.table))
cat(sprintf('%d of %d datasets are done with Nanopore of PacBio\n', sum(keep), length(keep)))

run.info.table..lr = run.info.table[keep,,drop = FALSE]

df..lr = cbind(run.info.table..lr[,c('bases'), drop = FALSE], as.data.frame(rep(1, nrow(run.info.table..lr))))
colnames(df..lr) = c('bases', 'runs_count')
df.aggr..lr = aggregate.data.frame(x = df..lr, by = list(BioProject.and.LibrarySource = sprintf('%s---%s', run.info.table..lr$BioProject, run.info.table..lr$LibrarySource)),  FUN = function(x){ sum (x) })

df.aggr.part..lr = t(as.data.frame(strsplit(df.aggr..lr$BioProject.and.LibrarySource, '---')))
colnames(df.aggr.part..lr) = c('BioProject', 'LibrarySource')

df.aggr..lr = cbind(df.aggr..lr, df.aggr.part..lr)
rownames(df.aggr..lr) = df.aggr..lr$BioProject.and.LibrarySource
colnames(df.aggr..lr) = c('BioProject.and.LibrarySource', 'bases', 'runs_count', 'BioProject', 'LibrarySource')

df.aggr.ext..lr = df.aggr..lr

for(field in c('ReleaseDate', 'CenterName')){
  all.data = run.info.table..lr[,field]
  names(all.data) = run.info.table..lr$BioProject
  all.data.nr = all.data[!duplicated(names(all.data))]
  all.data.nr.df = as.data.frame(all.data.nr[df.aggr..lr$BioProject])
  rownames(all.data.nr.df) = df.aggr..lr$BioProject.and.LibrarySource
  df.aggr.ext..lr = cbind(df.aggr.ext..lr, all.data.nr.df)
  colnames(df.aggr.ext..lr)[ncol(df.aggr.ext..lr)] = field
}

sel = which(sapply(df.aggr.ext..lr$CenterName, function(x) nchar(x) > 8))
df.aggr.ext..lr$CenterName[sel] = str_to_title(df.aggr.ext..lr$CenterName)[sel]




temp.aggr = aggregate.data.frame(x = run.info.table..lr[,c('Platform'), drop = FALSE],
                                 by = list(BioProject.and.LibrarySource = sprintf('%s---%s', run.info.table..lr$BioProject, run.info.table..lr$LibrarySource)),
                                 FUN = function(x){ 
                                   c_text = paste(levels(as.factor(x)), collapse = '; ')
                                   c_text = gsub('PACBIO_SMRT', 'PacBio', c_text, fixed = TRUE)
                                   c_text = gsub('OXFORD_NANOPORE', 'Nanopore', c_text, fixed = TRUE)
                                   c_text = gsub('ILLUMINA', 'Illumina', c_text, fixed = TRUE)
                                   c_text = gsub('BGISEQ', 'BGI', c_text, fixed = TRUE)
                                   c_text = gsub('ION_TORRENT', 'IonTorrent', c_text, fixed = TRUE)
                                })
rownames(temp.aggr) = temp.aggr$BioProject.and.LibrarySource
df.aggr.ext..lr = cbind(df.aggr.ext..lr, as.data.frame(temp.aggr[rownames(df.aggr.ext..lr), 'Platform']))
colnames(df.aggr.ext..lr)[length(df.aggr.ext..lr)] = 'Platforms'



temp.aggr = aggregate.data.frame(x = run.info.table..lr[,c('ScientificName'), drop = FALSE],
                                 by = list(BioProject.and.LibrarySource = sprintf('%s---%s', run.info.table..lr$BioProject, run.info.table..lr$LibrarySource)),
                                 FUN = function(x){ length(levels(as.factor(x))) })
rownames(temp.aggr) = temp.aggr$BioProject.and.LibrarySource
df.aggr.ext..lr = cbind(df.aggr.ext..lr, as.data.frame(temp.aggr[rownames(df.aggr.ext..lr), 'ScientificName']))
colnames(df.aggr.ext..lr)[length(df.aggr.ext..lr)] = 'N taxons'


temp.aggr = aggregate.data.frame(x = run.info.table..lr[,c('ScientificName'), drop = FALSE],
                                 by = list(BioProject.and.LibrarySource = sprintf('%s---%s', run.info.table..lr$BioProject, run.info.table..lr$LibrarySource)),
                                 FUN = function(c_taxons){ 
                                   c_table = table(c_taxons)
                                   c_table = c_table[rev(order(c_table))]
                                   c_top.taxons.info = sprintf('%s (%d)', names(c_table), c_table)
                                   
                                   if(length(c_table) > 3){
                                     paste(c(c_top.taxons.info[1:3],'....'), collapse = '; ')
                                   } else {
                                     paste(c_top.taxons.info, collapse = '; ')
                                   }
                                 })


rownames(temp.aggr) = temp.aggr$BioProject.and.LibrarySource
df.aggr.ext..lr = cbind(df.aggr.ext..lr, as.data.frame(temp.aggr[rownames(df.aggr.ext..lr), 'ScientificName']))
colnames(df.aggr.ext..lr)[length(df.aggr.ext..lr)] = 'Taxons'
df.aggr.ext..lr$ReleaseDate = as.Date(df.aggr.ext..lr$ReleaseDate)

omit = df.aggr.ext..lr$BioProject %in% ""
omit = omit | (df.aggr.ext..lr$bases == 0)
df.aggr.ext..lr = df.aggr.ext..lr[!omit,,drop = FALSE]
df.aggr.ext..lr = df.aggr.ext..lr[rev(order(df.aggr.ext..lr$bases)),]

write.table(df.aggr.ext..lr, sprintf('BioProject-centric.info.%s.tsv', descr), sep = '\t', quote = FALSE)



