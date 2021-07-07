import pandas as pd
import os
import re

configfile: "config/config.yml"

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
ref_path = os.path.join(
    config['ref']['root'], 
    "gencode-release-" + gencode,
    build,
    "dna"
    )
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = build + "." + assembly + ".fa"
ref_fagz = build + "." + assembly + ".fa.gz"

##########################################################
## Define all the required outputs from the setup steps ##
##########################################################

## Directories
raw_path = "data/raw/fastq"
trim_path = "data/trimmed/fastq"

FAGZ = [os.path.join(ref_path, ref_fagz)]
BOWTIEIDX = expand(
    [ref_path + "/bt2/{prefix}.{sub}.bt2"],
    prefix = build + "." + assembly,
    sub = ['1', '2', '3', '4', 'rev.1', 'rev.2']
)
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

ALL_OUTPUTS = []
ALL_OUTPUTS.extend([BOWTIEIDX])
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

#####################
## HiC-Pro outputs ##
#####################
hic_path = "data/hic"

# Required annotations during setup for HiC-Pro
chr_sizes = os.path.join("config", build + ".chr_sizes.tsv")
rs_frags = os.path.join(
    "config",
    build + "_" + config['hicpro']['enzyme'] + "_fragment.bed"
    )
REFS = [chr_sizes, rs_frags]

# Update the config file
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"
digest_script = "scripts/digest_genome.py"
# PROC_PAIRS = expand([hic_path + "/hic_results/data/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"],
#                     sample = samples)
# HIC_QC = expand([hic_path + "/hic_results/pic/{sample}"], sample = samples)
# VALID_PAIRS = expand([hic_path + "/hic_results/data/{sample}/{sample}_allValidPairs"],
#                        sample = samples)
# HIC_MAT = expand([hic_path + "/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
#                   sample = samples, bin = bins)
# HIC_BED = expand([hic_path + "/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_abs.bed"],
#                   sample = samples, bin = bins)
# MERGED_INT = expand([hic_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}{suffix}"],
#                     bin = bins, suffix = ['.matrix', '_abs.bed'])

ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend([hicpro_config, digest_script])
# ALL_OUTPUTS.extend(PROC_PAIRS)
# ALL_OUTPUTS.extend(HIC_QC)
# ALL_OUTPUTS.extend(VALID_PAIRS)
# ALL_OUTPUTS.extend(HIC_MAT)
# ALL_OUTPUTS.extend(HIC_BED)
# ALL_OUTPUTS.extend(MERGED_INT)


#####################
## Max HiC Outputs ##
#####################
# MAXHIC_INTERACTIONS = expand(["output/MaxHiC/merged/{bin}/{type}_interactions.txt.gz"],
#                              bin = bins, type = ['cis', 'trans'])
# ALL_OUTPUTS.extend(MAXHIC_INTERACTIONS)

#################################
## Detection of SuperEnhancers ##
#################################
# Additional Data for super-enhancers shold be in data/external/H3K27AC/bam
# h3k27ac_dir = 'data/external/H3K27AC'
# enh = ['T47D_H3K27Ac_E2', 'T47D_H3K27Ac_E2_DHT']
# SE_OUT =  expand(["output/HOMER/{sample}/{file}.tsv"],
#                  sample = enh, file = ['enhancers', 'superEnhancers'])
# ALL_OUTPUTS.extend(SE_OUT)


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
# include: "rules/merge_matrices.smk"
# include: "rules/maxhic.smk"

