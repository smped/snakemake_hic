library(readr)
library(dplyr)

args <- commandArgs(TRUE)
## Accept args in the following order
## 1 - The current bin size
## 2 - The input path
## 3 - The output path
## 4+ - All samples
## The last arg can thus be of variable length with the number of samples
## able to be calculatated from these arguments
# args <- c(
#     "5000",
#     "data/hic/hic_results/matrix/",
#     "data/hic/hic_results/matrix/merged/raw/",
#     "Veh1_S1",
#     "Veh2_S2",
#     "Veh4_S3"
# )
message("Supplied arguments are:\n", paste(args, collapse = "\n"))

## Construct the list of input matrices
bin <- args[[1]]
in_path <- args[[2]]
out_path <- file.path(args[[3]], bin)
stopifnot(length(args) > 3)
samples <- args[seq(4, length(args))]
mat_path <- file.path(
    in_path, samples, "raw", bin,
    paste0(samples, "_", bin, ".matrix")
)
stopifnot(all(file.exists(mat_path)))
## And the output file
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
mat_out <- file.path(
    out_path, paste0("merged_", bin, ".matrix")
)

## Construct the list of input bed files
bed_path <- file.path(
    in_path, samples, "raw", bin,
    paste0(samples, "_", bin, "_abs.bed")
)
stopifnot(all(file.exists(bed_path)))

## Check all bed files are identical
md5 <- vapply(bed_path, tools::md5sum, character(1))
data.frame(md5 = md5)
stopifnot(length(unique(md5)) == 1)

## Copy the bed file
bed_out <- file.path(
    out_path, paste0("merged_", bin, "_abs.bed")
)
file.copy(bed_path[[1]], bed_out, overwrite = TRUE)

## The matrix files may be >100 million lines long
## Try being lazy, but this may need to be chunked, or else we can try vroom
df <- lapply(
    mat_path,
    read_tsv,
    col_names = c("bin1", "bin2", "count"),
    col_types = "iii"
)
df <- bind_rows(df)
df <- arrange(df, bin1, bin2)
df <- group_by(df, bin1, bin2) # This blows the RAM out. Chunk from here
df <- summarise(df, count = sum(count), .groups = "drop")

## Export
write_tsv(df, mat_out, col_names = FALSE)

