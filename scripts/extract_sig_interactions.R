library(tidyverse)
library(magrittr)
library(yaml)
library(rtracklayer)
library(InteractionSet)
library(vroom)
library(glue)
library(pander)

args <- commandArgs(TRUE)
## This must simply be the bin size
stopifnot(length(args) == 1)
sz <- args[[1]]

## The FDR. Move to config later?
alpha <- 0.05

## Load the config file
config <- here::here("config/config.yml") %>%
  read_yaml()

## Define the correct seqinfo object to place into every InteractionSet
sq <- here::here(
  "output", config$ref$build,
  glue("{config$ref$build}.chr_sizes.tsv")
  ) %>% 
  read_tsv(col_names = c("seqnames", "length")) %>%
  with(
    Seqinfo(
      seqnames = seqnames,
      seqlengths = length,
      isCircular = rep(FALSE, length(seqnames)),
      genome = config$ref$build
    )
  )

## Use the bed file to define genomic bins
bins <- import.bed(
  here::here(
    "output/hic_pro/matrix",
    glue("merged_{sz}_abs.bed.gz")
  ),
  seqinfo = sq
)
## Load the cis interactions & filter based on the FDR
cis_int <- vroom(
  here::here(
    "output", "MaxHiC", sz, "cis_interactions.txt.gz"
  ),
  col_types = "dddd-d--"
)
fdr <- p.adjust(exp(-cis_int$neg_ln_p_val), "fdr")
keepInt <- fdr < alpha

## Setup the GInteractions & export the rds
gi <- GInteractions(
  anchor1 = bins[as.integer(cis_int$bin1ID[keepInt])],
  anchor2 = bins[as.integer(cis_int$bin2ID[keepInt])],
  counts = as.integer(cis_int$observed_interactions[keepInt]),
  exp_interactions = cis_int$exp_interactions[keepInt],
  p = exp(-cis_int$neg_ln_p_val)[keepInt],
  fdr = fdr[keepInt]
)
gi$bin_size <- Rle(sz, sum(keepInt))
gi$distance <- pairdist(gi)
gi$logObsExp <- log2(gi$counts / gi$exp_interactions)
seqinfo(gi) <- sq
write_rds(
  gi,
  here::here("output/MaxHiC", glue("gi_{sz}_cis.rds")),
  compress = "gz"
)

## Cleanup for times when the data is *big*
rm(cis_int)
rm(gi)
gc()

## Repeat for the trans interactions
trans_int <- vroom(
  here::here(
    "output", "MaxHiC", sz, "trans_interactions.txt.gz"
  ),
  col_types = "dddd-d--"
)
fdr <- p.adjust(exp(-trans_int$neg_ln_p_val), "fdr")
keepInt <- fdr < alpha
gi <- GInteractions(
  anchor1 = bins[as.integer(trans_int$bin1ID[keepInt])],
  anchor2 = bins[as.integer(trans_int$bin2ID[keepInt])],
  counts = as.integer(trans_int$read_count[keepInt]),
  p = exp(-trans_int$neg_ln_p_val)[keepInt],
  fdr = fdr[keepInt]
)
gi$bin_size <- Rle(sz, sum(keepInt))
seqinfo(gi) <- sq
write_rds(
  gi,
  here::here("output/MaxHiC", glue("gi_{sz}_trans.rds")),
  compress = "gz"
)

