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

ALL_OUTPUTS = []

## Genome & index
FAGZ = [os.path.join(ref_path, ref_fagz)]
ALL_OUTPUTS.extend(FAGZ)

## Directories for data
raw_path = "data/raw/fastq"
trim_path = "data/trimmed/fastq"
FQC_OUTS = expand(
    ["data/{step}/FastQC/{sample}_{reads}_fastqc.{suffix}"],
    suffix = ['zip', 'html'],
    reads = ['R1', 'R2'],
    sample = list(df['path']),
    step = ['raw', 'trimmed']
)
TRIM_OUTS = expand(
    [trim_path + "/{sample}_{reads}{suffix}"],
    sample = list(df['path']),
    suffix = suffix,
    reads = ['R1', 'R2']
)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

#####################
## HiC-Pro outputs ##
#####################
hic_data_path = "data/hic"
bowtie_data_path = os.path.join(hic_data_path, "bowtie_results")

# Required annotations during setup for HiC-Pro
chr_sizes = os.path.abspath(
    os.path.join(
        "output", 
        build,
        build + ".chr_sizes.tsv"
    )
)
rs_frags = os.path.abspath(
    os.path.join(
        "output", 
        build,
        build + "_" + config['hicpro']['enzyme'] + "_fragment.bed"
        )
)
ALL_OUTPUTS.extend([rs_frags + ".gz"])

# Update the config file
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"
ALL_OUTPUTS.extend([hicpro_config])

# Get the mappings & QC
HIC_QC = expand(
    [hic_data_path + "/hic_results/pic/{sample}/plotMapping_{sample}.pdf"],
    sample = samples
    )
HIC_MERGE_PAIRS = expand(
    [hic_data_path + "/hic_results/data/{sample}/{sample}.allValidPairs"],
    sample = samples
    )
HIC_MERGE_STAT = expand(
    [hic_data_path + "/hic_results/stats/{sample}/{sample}{suffix}"],
    sample = samples,
    suffix = ['.mRSstat', '.mpairstat', read_ext[0] + ".mmapstat", read_ext[1] + ".mmapstat", "_allValidPairs.mergestat"]
    )
HIC_CONTACT_MAPS = expand(
    [hic_data_path + "/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}{suffix}"],
    sample = samples,
    bin = bins,
    suffix = ['.matrix', '_abs.bed']
    )
HIC_CONTACT_PICS = expand(
    [hic_data_path + "/hic_results/pic/{sample}/plot{file}_{sample}.pdf"],
    sample = samples,
    file = ['HiCContactRanges', 'HiCFragmentSize', 'HiCFragment', 'MappingPairing']
    )
ALL_OUTPUTS.extend([HIC_QC])
ALL_OUTPUTS.extend([HIC_MERGE_PAIRS, HIC_MERGE_STAT])
ALL_OUTPUTS.extend([HIC_CONTACT_MAPS])

## Merge the interaction matrices
MERGED_INT = expand(
    [hic_data_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}{suffix}"],
    bin = bins,
    suffix = ['.matrix', '_abs.bed']
    )
ALL_OUTPUTS.extend(MERGED_INT)


#####################
## Max HiC Outputs ##
#####################
MAXHIC_INTERACTIONS = expand(
    ["output/MaxHiC/merged/{bin}/{type}_interactions.txt.gz"],
    bin = bins, type = ['cis', 'trans']
    )
ALL_OUTPUTS.extend(MAXHIC_INTERACTIONS)


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
