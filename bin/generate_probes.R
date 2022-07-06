#!/usr/bin/env Rscript

library(tidyverse)
library(utils)
require(scales)
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

print(n.name)
print(v.name)

# # Read in the tables with the tandem repeat units:

# tr.units.n = read.delim("N_top10pc_tab_bycol.tsv", header = F, sep = "\t", stringsAsFactors = F) 
# names(tr.units.n) = c("contig", "length", "mean.coverage", "repr", "unit.seq")

# tr.units.v = read.delim("V_top10pc_tab_bycol.tsv", header = F, sep = "\t", stringsAsFactors = F)
# names(tr.units.v) = c("contig", "length", "mean.coverage", "repr", "unit.seq")

# tr.units.n = tr.units.n %>%
#   dplyr::arrange(unit.seq) %>%
#   dplyr::filter(nchar(unit.seq) > 6) %>%
#   group_by(unit.seq) %>%
#   do(dplyr::arrange(., desc(.$mean.coverage))) %>%
#   do(head(., 1)) %>%
#   ungroup()

# tr.units.v = tr.units.v %>%
#   dplyr::arrange(unit.seq) %>%
#   dplyr::filter(nchar(unit.seq) > 6) %>%
#   group_by(unit.seq) %>%
#   do(dplyr::arrange(., desc(.$mean.coverage))) %>%
#   do(head(., 1)) %>%
#   ungroup()

# # Calculate the Levenshteine distance between all K-bp substrings from the repeat units of D. nairensis and D. valentini.

# probe.length = 30

# # Generate all k-mers from the string s; if k is greater than the length of s then replicate s.
# take_all_substrings = function(s, k) {
#   if (nchar(s) < k) {
#     s = paste0(rep(s, k %/% nchar(s) + 1), collapse = '')
#   }
#   return(unlist(purrr::map(1:(nchar(s) - k + 1), 
#                              function(i) stringr::str_sub(s, i, i + k - 1))))
# }

# # Make a table of all probes of length k from a given repeat sequence found in a given contig.
# generate_substring_df = function(x, k) {
#   s = x[1, "unit.seq"]
#   df = data.frame(probe.seq = take_all_substrings(s, k),
#                   stringsAsFactors = FALSE)
#   df$contig = x$contig
#   df$length = x$length
#   df$mean.coverage = x$mean.coverage
#   df$repr = x$repr
#   return(df)
# }

# # Check if more than a half of the length of a probe is not a microsatellite array sequence.
# # Length of a probe is assumed to be 30 bp.
# rm_microsatellites = function(probe.seq) {
#   for(unit.length in 2:6) {
#     unit.num = 6 - unit.length + 2 # The number of repeat units that occupy about a half of the 30-bp probe.
#     coord.last = probe.length - unit.num * unit.length + 1
#     delta.seq = seq(0, unit.length * (unit.num - 1), by = unit.length)
#     for(coord in 1:coord.last) {
#       unique.units = unique(
#         unlist(
#         lapply(delta.seq, function(delta) {
#         stringr::str_sub(probe.seq, coord + delta, coord + delta + unit.length - 1)
#       })))
#       if (length(unique.units) == 1) {
#         return(FALSE)
#       }
#     }
#   }
#   return(TRUE)
# }

# # Generate probes.

# tr.probes.n = tr.units.n %>%
#   group_by(contig, unit.seq) %>%
#   do(generate_substring_df(., probe.length)) %>%
#   ungroup() %>%
#   distinct() %>%
#   group_by(probe.seq, contig) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "AAAAAA")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "TTTTTT")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "GGGGGG")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "CCCCCC")) %>%
#   dplyr::filter(rm_microsatellites(probe.seq)) %>%
#   ungroup()

# tr.probes.v = tr.units.v %>%
#   group_by(contig, unit.seq) %>%
#   do(generate_substring_df(., probe.length)) %>%
#   ungroup() %>%
#   distinct() %>%
#   group_by(probe.seq, contig) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "AAAAAA")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "TTTTTT")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "GGGGGG")) %>%
#   dplyr::filter(!stringr::str_detect(probe.seq, "CCCCCC")) %>%
#   dplyr::filter(rm_microsatellites(probe.seq)) %>%
#   ungroup()

# # Calculate Levenshtein distance between each pair of D. raddei nairensis and D. valentini probes. 
# names(tr.probes.n) = c("unit.seq.n", "probe.seq.n", "contig.n", "length.n", "mean.coverage.n", "repr.n")
# names(tr.probes.v) = c("unit.seq.v", "probe.seq.v","contig.v", "length.v", "mean.coverage.v", "repr.v")
# tr.probes.n$dummy.col = 1
# tr.probes.v$dummy.col = 1
# tr.probes.all = tr.probes.n %>%
#   inner_join(tr.probes.v, by = ("dummy.col" = "dummy.col")) %>%
#   group_by(unit.seq.n, contig.n, probe.seq.n, unit.seq.v, contig.v, probe.seq.v) %>%
#   dplyr::mutate(edit.dist = utils::adist(probe.seq.n, probe.seq.v)) %>%
#   ungroup() %>%
#   dplyr::select(-dummy.col)

# # Calculate the minimal distance between each probe for one species and the whole set of probes for the other species.
# tr.probes.min.all = tr.probes.all %>%
#   group_by(probe.seq.n) %>%
#   mutate(min.dist.n = min(edit.dist)) %>%
#   ungroup() %>%
#   group_by(probe.seq.v) %>%
#   mutate(min.dist.v = min(edit.dist)) %>%
#   ungroup() %>%
#   group_by(probe.seq.n) %>%
#   mutate(min.n = (edit.dist == min.dist.n)) %>%
#   ungroup() %>%
#   group_by(probe.seq.v) %>%
#   mutate(min.v = (edit.dist == min.dist.v)) %>%
#   ungroup()

# # For each probe for D. raddei nairensis select the closest probe for D. valentini.
# tr.probes.n.arranged = tr.probes.min.all %>%
#   filter(min.n) %>%
#   dplyr::arrange(desc(edit.dist)) %>%
#   dplyr::filter(min.dist.n > 10)

# # For each probe for D. valentini select the closest probe for D. raddei nairensis.
# tr.probes.v.arranged = tr.probes.min.all %>%
#   filter(min.v) %>%
#   dplyr::arrange(desc(edit.dist)) %>%
#   dplyr::filter(min.dist.v > 10)

# write.table(tr.probes.n.arranged,
#             "N_tr_probes_by_distance_geq10.tsv",
#             quote = F,
#             sep = "\t",
#             row.names = F)

# write.table(tr.probes.v.arranged,
#             "V_tr_probes_by_distance_geq10.tsv",
#             quote = F,
#             sep = "\t",
#             row.names = F)
