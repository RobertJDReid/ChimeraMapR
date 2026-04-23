library(data.table)
library(ggplot2)

d1 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0500_loess_fits_2026-04-20.csv")
d2 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0250_loess_fits_2026-04-20.csv")
d3 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0100_loess_fits_2026-04-20.csv")
d4 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0050_loess_fits_2026-04-20.csv")
d5 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0025_loess_fits_2026-04-20.csv")
d6 <- fread("~/Dropbox/__Work/__Lab-Notebook/2026/04/2026-04-16 Read analysis downsample/loess/Sample_0010_loess_fits_2026-04-20.csv")





all_fits <- rbindlist(list(d1,d2,d3,d4,d5,d6), fill = TRUE)

plot_chrom = "S288C_chrXVI"

ggplot(all_fits[chrom == plot_chrom],
       aes(x = pos_kb, y = uniform_fit,
           fill = factor(sample_name,levels = 
                           c("Sample_0500","Sample_0250",
                             "Sample_0100","Sample_0050",
                             "Sample_0025","Sample_0010")))) +
  geom_ribbon(aes(ymin=0,ymax=uniform_fit)) +
  geom_line(linewidth = 0.1,color="black") +
  theme_bw() +
  scale_fill_viridis_d(option = "turbo",begin = 0.85,end = 0.15) +
  #scale_fill_grey(start = 0.9, end = 0) +
  xlab(paste0(plot_chrom," Position (Kb)")) +
  ylab("LOESS fit") +
  theme(
    legend.title = element_blank()
  )
