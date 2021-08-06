import pandas as pd
import os
import re
import subprocess

configfile: "config/config.yml"

#####################################
## Check that HiC-Pro is installed ##
#####################################
# Get the path if installed
hic_check = subprocess.run(
    ['which', 'HiC-Pro'],
    universal_newlines=True,
    check=True,
    stdout=subprocess.PIPE
)
hic_path = hic_check.stdout
rc = hic_check.returncode
if not rc == 0:
    print("HiC-Pro not installed")
    sys.exit(1)

##################
# Define Samples #
##################
df = pd.read_table(config["samples"])
samples = list(set(df['sample']))
cols = list(df.columns)
df['path'] = df[cols].apply(lambda row: '/'.join(row.values.astype(str)), axis=1)
suffix = config['suffix']
# Can this be done better?
read_ext = [config['hicpro']['pair1_ext'], config['hicpro']['pair2_ext']]

#################################
## Variables for the reference ##
#################################

build = config['ref']['build']
genbank = config['ref']['genbank']
gencode = str(config['ref']['gencode'])
ref_path = os.path.abspath(
    os.path.join(
        config['ref']['root'],
        "gencode-release-" + gencode,
        build,
        "dna"
    )
)
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = build + "." + assembly + ".fa"
ref_fagz = build + "." + assembly + ".fa.gz"

##########################################################
## Define all the required outputs from the setup steps ##
##########################################################

## Directories for data
raw_path = "data/raw/fastq"
trim_path = "data/trimmed/fastq"
ALL_OUTPUTS = []

#####################
## HiC-Pro outputs ##
#####################
hic_data_path = "data/hic"
hic_output_path = os.path.join("output", "hic_pro")
bowtie_data_path = os.path.join(hic_data_path, "bowtie_results")

# Required annotations during setup for HiC-Pro
chr_sizes = os.path.abspath(
    os.path.join("output", build, build + ".chr_sizes.tsv")
)
rs_frags = os.path.abspath(
    os.path.join(
        "output",
        build,
        build + "_" + config['hicpro']['enzyme'] + "_fragment.bed"
        )
)
fragment_lengths = os.path.join("output", build, "fragment_length.counts")
## We can exclude them from the required output as they're only
## needed as input for other rules.
ALL_OUTPUTS.extend([rs_frags + ".gz"])

# Update the config file
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"

## Merge the interaction matrices
MERGED_INT = expand(
    [
        os.path.join(
            hic_output_path, "matrix", "merged_{bin}{suffix}.gz"
        )
    ],
    bin = bins,
    suffix = ['.matrix', '_abs.bed']
    )
ALL_OUTPUTS.extend(MERGED_INT)

## Collect the final required outputs
HIC_VALID_PAIRS = expand(
    [
        os.path.join(
            hic_output_path, "allValidPairs", "{sample}.allValidPairs.gz"
        )
    ],
    sample = samples
    )
HIC_MATRICES = expand(
    [
        os.path.join(
            hic_output_path, "matrix", "raw", "{bin}",
            "{sample}_{bin}.matrix.gz"
        )
    ],
    sample = samples,
    bin = bins
)
HIC_STATS = expand(
    [os.path.join(hic_output_path, "stats", "{sample}", "{sample}{suffix}")],
    sample = samples,
    suffix = [
        '.mRSstat', read_ext[0] + ".mmapstat", read_ext[1] + ".mmapstat",
        "_allValidPairs.mergestat", ".mpairstat"
    ]
)
HIC_QC_PDF = expand(
    [
        os.path.join(
            hic_data_path, "hic_results", "pic", "{sample}",
            "plot{file}_{sample}.pdf"
        )
    ],
    sample = samples,
    file = ['HiCFragment', 'MappingPairing', 'Mapping']
)
ALL_OUTPUTS.extend(
    [HIC_VALID_PAIRS, HIC_MATRICES, HIC_STATS, HIC_QC_PDF]
)


#####################
## Max HiC Outputs ##
#####################
# MAXHIC_INTERACTIONS = expand(
#     ["output/MaxHiC/{bin}/{type}_interactions.txt.gz"],
#     bin = bins, type = ['cis', 'trans']
#     )
# ALL_OUTPUTS.extend(MAXHIC_INTERACTIONS)


##########################
## QC Workflowr Reports ##
##########################
wd = os.getcwd()
rproj = os.path.basename(wd) + ".Rproj"
# WFLOW_OUT = expand(
#     ["docs/{file}.html"],
#     file = ['index', 'qc_raw', 'qc_trimmed', 'qc_hic', 'define_interactions']
# )
# ALL_OUTPUTS.extend([rproj])
# ALL_OUTPUTS.extend(WFLOW_OUT)

#####################
## Rules & Outputs ##
#####################

rule all:
    input:
        ALL_OUTPUTS

include: "rules/reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/hicpro.smk"
include: "rules/merge_matrices.smk"
include: "rules/maxhic.smk"
include: "rules/collect_output.smk"
include: "rules/workflowr.smk"
