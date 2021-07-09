library(yaml)
library(stringr)
library(magrittr)

## These cannot be named so need to be carefully specified in order
## 1 - The Bowtie2 Index path
## 2 - The chr_sizes file
## 3 - The genome fragments file
## 4 - The path to the HiC-Pro template
## 5 - The location of the output file
args <- commandArgs(TRUE)
stopifnot(length(args) == 5)

## Load the local YAML config file
config <- read_yaml(here::here("config/config.yml"))

## Import the HiC-Pro config file from the installed version
template <- args[[4]]
stopifnot(file.exists(template))
orig <- readLines(template)

## Modify & write the file
orig %>%
  str_replace("^N_CPU.+", paste("N_CPU =", config$hicpro$ncpu)) %>%
  str_replace("^SORT_RAM.+", paste("SORT_RAM =", config$hicpro$sort_ram)) %>%
  str_replace("^PAIR1_EXT.+", paste("PAIR1_EXT =", config$hicpro$pair1_ext)) %>%
  str_replace("^PAIR2_EXT.+", paste("PAIR2_EXT =", config$hicpro$pair2_ext)) %>%
  str_replace("^FORMAT.+", paste0("FORMAT = phred", config$hicpro$phred)) %>%
  str_replace("^MIN_MAPQ.+", paste("MIN_MAPQ =", config$hicpro$min_mapq)) %>%
  str_replace("^BOWTIE2_IDX_PATH.+", paste("BOWTIE2_IDX_PATH =", args[[1]])) %>%
  str_replace(
    "^REFERENCE_GENOME.+",
    paste0("REFERENCE_GENOME = ", config$ref$build, ".", config$ref$assembly, ".genome")
  ) %>%
  str_replace("^GENOME_SIZE.+", paste("GENOME_SIZE =", args[[2]])) %>%
  str_replace("^GENOME_FRAGMENT.+", paste("GENOME_FRAGMENT =", args[[3]])) %>%
  str_replace("^LIGATION_SITE.+", paste("LIGATION_SITE =", config$hicpro$ligation_site)) %>%
  str_replace("^BIN_SIZE.+", paste("BIN_SIZE =", config$hicpro$bin_size)) %>%
  str_replace("^MATRIX_FORMAT.+", paste("MATRIX_FORMAT =", config$hicpro$matrix_format)) %>%
  writeLines(con = args[[5]])

