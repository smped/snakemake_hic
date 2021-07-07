# snakemake_hic

This is a template repository for running HiC-Pro under snakemake.
The steps currently implemented are:

1. Download and index the genome
2. Trim raw data using `AdapterRemoval`
3. Run `FastQC` on raw and trimmed data
4. Prepare Files for HiC-Pro
    + Organise genome annotations 
    + Define restriction fragments 
    + Update `hicpro-config.txt` based on `config.yml`

## Setup

HiC-Pro must be installed and visible to this workflow.
Please ensure you have a working installation on your system.
This workflow is built around v3.0.0

- Data must be placed in `data/raw/fastq`
    + Each sample/replicate needs to go in a separate folder
- Edit `config/samples.tsv` to reflect the sample names and file names. 
    + Data is assumed to be paired, so only one file needs to be listed per directory
- Edit `config/config.yml`

## Testing

- Once all edits are performed and data is uploaded, run `snamekmake -n` as a dry run to check everything will work as expected
- If the dry run is successful, create the conda environments `snakemake --use-conda --conda-create-envs-only`. This may take a while.

Create the rulegraph to check everything is as expected

```
snakemake --rulegraph > .snakemake/rulegraph.dot
dot -Tpdf .snakemake/rulegraph.dot > .snakemake/rulegraph.pdf
```

## Execution

```
snakemake \
    --use-conda \
    --notemp \
    --cores 12
```

