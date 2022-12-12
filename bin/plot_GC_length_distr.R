#!/usr/bin/env Rscript

library(tidyverse)
require(scales)
stringsAsFactors = F

# Upload tables of contigs and repeat units:

contigs.n <- read.delim('GC_R/contigs_N_tab.tsv', header = F,  sep="\t", stringsAsFactors = F)
names(contigs.n) = c("contig", "length", "mean.coverage", "repr", "unit.seq")
contigs.n <- contigs.n %>% drop_na()

contigs.v <- read.delim('GC_R/contigs_V_tab.tsv', header = F,  sep="\t", stringsAsFactors = F) 
names(contigs.v) = c("contig", "length", "mean.coverage", "repr", "unit.seq")
contigs.v <- contigs.v %>% drop_na()

all.repeats.n <- read.delim('GC_R/N_all_tab_bycol.tsv', header = F,  sep="\t", stringsAsFactors = F) 
names(all.repeats.n) = c("contig", "length", "mean.coverage", "repr", "unit.seq")

all.repeats.v <- read.delim('GC_R/V_all_tab_bycol.tsv', header = F,  sep="\t", stringsAsFactors = F)
names(all.repeats.v) = c("contig", "length", "mean.coverage", "repr", "unit.seq")
 
# Define a function to calcualte GC-content:

get_GC_content <- function(unit_seq) {
  num_g <- str_count(unit_seq, "G")
  num_c <- str_count(unit_seq, "C")
  gc_content <- (num_g + num_c) / str_length(unit_seq) * 100 
  return(as.numeric(gc_content))
}

# Calculate the GC-content and length of contigs and repeat units:

contigs.n <- contigs.n %>%
  mutate(GC = get_GC_content(unit.seq)) %>%
  mutate(unit.length = nchar(unit.seq)) 
contigs.v <- contigs.v %>%
  mutate(GC = get_GC_content(unit.seq)) %>%
  mutate(unit.length = nchar(unit.seq)) 

all.repeats.n <- all.repeats.n %>%
  mutate(GC = get_GC_content(unit.seq)) %>%
  mutate(unit.length = nchar(unit.seq)) 
all.repeats.v <- all.repeats.v %>%
  mutate(GC = get_GC_content(unit.seq)) %>%
  mutate(unit.length = nchar(unit.seq))

# Create an extended table for contigs to normalize on the read coverage:

contigs.n.subset = contigs.n %>%
  dplyr::select(mean.coverage, GC, unit.length) 

contigs.n.extended = as.data.frame(lapply(contigs.n.subset, rep, contigs.n.subset$mean.coverage)) %>%
  mutate(species = "D. raddei nairensis")

contigs.v.subset = contigs.v %>%
  dplyr::select(mean.coverage, GC, unit.length) 

contigs.v.extended = as.data.frame(lapply(contigs.v.subset, rep, contigs.v.subset$mean.coverage)) %>%
  mutate(species = "D. valentini")

contigs.extended = dplyr::bind_rows(contigs.n.extended, contigs.v.extended) %>% 
  mutate(species = factor(species, levels = c("D. raddei nairensis", "D. valentini")))

# Create an extended table for repeat units to normalize on the read coverage:

all.repeats.n.subset = all.repeats.n %>%
  dplyr::select(mean.coverage, GC, unit.length) 

all.repeats.n.extended = as.data.frame(lapply(all.repeats.n.subset, rep, all.repeats.n.subset$mean.coverage)) %>%
  mutate(species = "D. raddei nairensis")

all.repeats.v.subset = all.repeats.v %>%
  dplyr::select(mean.coverage, GC, unit.length) 

all.repeats.v.extended = as.data.frame(lapply(all.repeats.v.subset, rep, all.repeats.v.subset$mean.coverage)) %>%
  mutate(species = "D. valentini")

all.repeats.extended = dplyr::bind_rows(all.repeats.n.extended, all.repeats.v.extended) %>% 
  mutate(species = factor(species, levels = c("D. raddei nairensis", "D. valentini")))

# Plot GC-content distribution of contigs:

plot.contigs.gc = ggplot() +
  geom_histogram(data = contigs.n.extended,
                 alpha = 1,
                 aes(x = GC,
                     y = stat(count / sum(count)),
                     fill = "D. raddei nairensis"),
                 bins = 100) +
  geom_histogram(data = contigs.v.extended,
                 alpha = 0.7,
                 aes(x = GC,
                     y = stat(count / sum(count)),
                     fill = "D. valentini"),
                 bins = 100) +
  scale_fill_manual(values = c("D. raddei nairensis" = "#da2c38",
                               "D. valentini" = "#00d103"),
                    name = "Species") +
  xlab("Contig GC-content") +
  ylab("Normalized sum coverage") +
  scale_x_continuous(breaks = seq(0, 100, 5),
                     limits = c(-5, 105)) +
  scale_y_continuous(breaks = seq(0, 0.08, 0.01),
                     limits = c(-0.001, 0.08)) +
  theme_classic() +
  theme(legend.position = c(.95, .70),
        legend.justification = c("right", "top"),
        text = element_text(size = 15))

ggsave("GC_R/contig_gc_content.pdf",
       plot = plot.contigs.gc,
       width = 15,
       height = 10,
       scale = 1, 
       units = c("cm"),
       limitsize = FALSE)

# Plot GC-content distribution of repeat units:

plot.repeats.gc = ggplot() +
  geom_histogram(data = all.repeats.n.extended,
                 alpha = 1,
                 aes(x = GC, y = stat(count / sum(count)),
                     fill = "D. raddei nairensis"),
                 binwidth = 5) +
  geom_histogram(data = all.repeats.v.extended,
                 alpha = 0.7,
                 aes(x = GC, y = stat(count / sum(count)),
                     fill = "D. valentini"),
                 binwidth = 5) +
  scale_fill_manual(values = c("D. raddei nairensis" = "#da2c38",
                               "D. valentini" = "#00d103"),
                    name = "Species") +
  xlab("Unit GC-content") +
  ylab("Normalized sum coverage") +
  scale_x_continuous(breaks = seq(0, 100, 5))+
  scale_y_continuous(breaks = seq(0, 0.3, 0.1),
                     limits = c(0, 0.3)) +
  theme_classic() +
  theme(legend.position = c(.37, .90),
        legend.justification = c("right", "top"),
        text = element_text(size = 15))

ggsave("GC_R/units_gc_content.pdf",
       plot = plot.repeats.gc,
       width = 15,
       height = 10,
       scale = 1,
       units = c("cm"),
       limitsize = FALSE)

# Plot sequence length distribution of contigs:

plot.contigs.length = ggplot(data = contigs.extended,
                             aes(x = species,
                                 y = unit.length,
                                 fill = species)) +
  geom_violin() +
  stat_summary(fun = median,
               geom = "point",
               size = 2,
               color = "black") + 
  scale_fill_manual(values = c("#da2c38",
                             alpha("#00d103", 0.7))) +
  ylab("Contig length") +
  scale_y_continuous(breaks = seq(0, 1100, 100),
                     limits = c(0, 1100)) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 15),
        axis.title.x = element_blank())

ggsave("GC_R/contigs_length.pdf",
       plot = plot.contigs.length,
       width = 15,
       height = 10,
       scale = 1,
       units = c("cm"),
       limitsize = FALSE)

# Plot sequence length distribution of repeat units:

plot.repeats.length = ggplot() +
  geom_histogram(data = all.repeats.n.extended,
                 alpha = 1,
                 aes(x = unit.length, y = stat(count / sum(count)),
                     fill = "D. raddei nairensis"),
                 binwidth = 10) +
  geom_histogram(data = all.repeats.v.extended,
                 alpha = 0.7,
                 aes(x = unit.length, y = stat(count / sum(count)),
                     fill = "D. valentini"),
                 binwidth = 10) +
  scale_fill_manual(values = c("D. raddei nairensis" = "#da2c38",
                               "D. valentini" = "#00d103"),
                    name = "Species") +
  xlab("Unit length") +
  ylab("Normalized sum coverage") +
  scale_x_continuous(breaks = seq(0, 180, 10)) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.1),
                     limits = c(0, 0.8)) +
  theme_classic() +
  theme(legend.position = c(.95, .70),
        legend.justification = c("right", "top"),
        text = element_text(size = 15))

ggsave("GC_R/units_length.pdf",
       plot = plot.repeats.length,
       width = 20,
       height = 10,
       scale = 1,
       units = c("cm"),
       limitsize = FALSE)

