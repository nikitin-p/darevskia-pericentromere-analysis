#!/usr/bin/env Rscript

# library("readxl")
library(tidyverse)
library(utils)
stringsAsFactors = F

# Parse input file names

args = commandArgs()
arg1 = args[8]
arg2 = args[9]

if (stringr::str_detect(arg1, fixed("N_"))) {
  n.name = arg1
  v.name = arg2
} else {
  n.name = arg2
  v.name = arg1
}

# Read in the tables with the tandem repeat units:

tr.units.n = read.delim(n.name, header = F, sep = "\t", stringsAsFactors = F)
names(tr.units.n) = c("contig", "length", "mean.coverage", "repr", "unit.seq")
tr.units.v = read.delim(v.name, header = F, sep = "\t", stringsAsFactors = F)
names(tr.units.v) = c("contig", "length", "mean.coverage", "repr", "unit.seq")

# Filter the repeat units:

tr.units.n = tr.units.n %>%
  dplyr::arrange(unit.seq) %>%
  dplyr::filter(nchar(unit.seq) > 6) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "AAAAAA")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "TTTTTT")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "GGGGGG")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "CCCCCC")) %>%
  group_by(unit.seq) %>%
  do(dplyr::arrange(., desc(.$mean.coverage))) %>%
  do(head(., 1)) %>%
  ungroup()

tr.units.v = tr.units.v %>%
  dplyr::arrange(unit.seq) %>%
  dplyr::filter(nchar(unit.seq) > 6) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "AAAAAA")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "TTTTTT")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "GGGGGG")) %>%
  dplyr::filter(!stringr::str_detect(unit.seq, "CCCCCC")) %>%
  group_by(unit.seq) %>%
  do(dplyr::arrange(., desc(.$mean.coverage))) %>%
  do(head(., 1)) %>%
  ungroup()

# Calculate the Levenshtine distance between each pair of repeat units from D. nairensis and from D. valentini:

names(tr.units.n) = c("contig.n", "length.n", "mean.coverage.n", "repr.n", "unit.seq.n")
names(tr.units.v) = c("contig.v", "length.v", "mean.coverage.v", "repr.v", "unit.seq.v")
tr.units.n$dummy.col = 1
tr.units.v$dummy.col = 1
tr.units.all = tr.units.n %>%
  inner_join(tr.units.v, by = ("dummy.col" = "dummy.col")) %>%
  group_by(contig.n, unit.seq.n, contig.v, unit.seq.v) %>%
  dplyr::mutate(edit.dist = utils::adist(unit.seq.n, unit.seq.v)) %>%
  ungroup() %>%
  group_by(contig.n) %>%
  dplyr::mutate(min.n = (edit.dist == min(edit.dist))) %>%
  ungroup() %>%
  group_by(contig.v) %>%
  dplyr::mutate(min.v = (edit.dist == min(edit.dist))) %>%
  ungroup() %>%
  dplyr::select(-dummy.col)

# Generate tables of repeat units, arranged by distance, for D. nairensis and for D. valentini.

# For D. nairensis:

tr.units.n.arranged = tr.units.all %>%
  filter(min.n) %>%
  dplyr::select(contig.n, mean.coverage.n, unit.seq.n, edit.dist) %>%
  dplyr::arrange(desc(edit.dist)) %>%
  distinct() 

write.table(tr.units.n.arranged,
            "N_tr_units_by_distance_geq10.tsv",
            quote = F,
            sep = "\t",
            row.names = F)


# For D. valentini:

tr.units.v.arranged = tr.units.all %>%
  filter(min.v) %>%
  dplyr::select(contig.v, mean.coverage.v, unit.seq.v, edit.dist) %>%
  dplyr::arrange(desc(edit.dist)) %>%
  distinct()

write.table(tr.units.v.arranged,
            "V_tr_units_by_distance_geq10.tsv",
            quote = F,
            sep = "\t",
            row.names = F)