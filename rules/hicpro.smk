rule find_rs_fragments:
    input: rules.unzip_reference.output
    output:
        rs = rs_frags
    params:
        enzyme = config['hicpro']['enzyme'],
        script = os.path.join(hic_path, "utils", "digest_genome.py")
    threads: 1
    conda: "../envs/hicpro.yml"
    shell:
        """
        # Run the python script
        python {params.script} \
          -r {params.enzyme} \
          -o {output.rs} \
          {input}
        """

rule make_hicpro_config:
    input:
        idx = rules.bowtie2_index.output[0],
        rs = rs_frags,
        chr_sizes = chr_sizes
    output:
        hicpro_config
    conda: "../envs/stringr.yml"
    threads: 1
    shell:
        """
        IDX=$(dirname {input.idx})
        Rscript --vanilla \
          scripts/write_hicpro_config.R \
          $IDX \
          {input.chr_sizes} \
          {input.rs} \
          {output}
        """

# rule hicpro_mapping:
#     input:
#         config = hicpro_config,
#         files = expand([trim_path + "/{{sample}}/{{sample}}{reads}{suffix}"],
#                        suffix = suffix, reads = read_ext)
#     output:
#         bam = temp(
#             expand(
#                 ["data/hic/bowtie_results/bwt2/{{sample}}/{{sample}}{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
#                 reads = read_ext
#             )
#         ),
#         glob = temp(
#             expand(
#                 [hic_data_path + "/bowtie_results/bwt2_global/{{sample}}/{{sample}}{reads}_" + build + "." + assembly + ".bwt2glob.{suffix}"],
#                 reads = read_ext,
#                 suffix = ['bam', 'unmap.fastq', 'unmap_trimmed.fastq']
#             )
#         ),
#         local = temp(
#             expand(
#                 [hic_data_path + "/bowtie_results/bwt2_local/{{sample}}/{{sample}}{reads}_" + build + "." + assembly + ".bwt2glob.unmap_bwt2loc.bam"],
#                 reads = read_ext
#              )
#         )
#     params:
#         indir = trim_path,
#         outdir = hic_data_path
#     log: "logs/hicpro/hicpro_mapping_{sample}.log"
#     threads: config['hicpro']['ncpu']
#     shell:
#         """
#         ######################################
#         ## Specific to phoenix for now only ##
#         ######################################
#         ## Load modules
#         module load HiC-Pro/2.9.0-foss-2016b

#         ## Remove any existing data as leaving this here causes HicPro to
#         ## make an interactive request. Piping `yes` into HicPro may be the
#         ## source of some current problems
#         if [[ -d {params.outdir} ]]; then
#           rm -rf {params.outdir}
#         fi

#         ## Run HiC-pro
#         HiC-Pro \
#           -s mapping \
#           -c {input.config} \
#           -i {params.indir} \
#           -o {params.outdir} &> {log}
#         """

# rule hicpro_proc:
#     input:
#         config = hicpro_config,
#         files = expand(
#                 [hic_data_path + "/bowtie_results/bwt2/{{sample}}/{{sample}}{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
#                 reads = read_ext
#             )
#     output:
#         bam = temp(hic_data_path + "/bowtie_results/bwt2/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.bam"),
#         pairs = hic_data_path + "/hic_results/data/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"
#     params:
#         indir = hic_data_path + "/bowtie_results/bwt2",
#         outdir = hic_data_path
#     log: "logs/hicpro/hicpro_proc_{sample}.log"
#     threads: config['hicpro']['ncpu']
#     shell:
#         """
#         ######################################
#         ## Specific to phoenix for now only ##
#         ######################################
#         ## Load modules
#         module load HiC-Pro/2.9.0-foss-2016b

#         ##Run HiC-pro responding to yes to any interactive requests
#         HiC-Pro \
#           -s proc_hic \
#           -c {input.config} \
#           -i {params.indir} \
#           -o {params.outdir} &> {log}
#         """

# rule hicpro_qc:
#     input:
#         config = hicpro_config,
#         files = expand(
#                 [hic_data_path + "/bowtie_results/bwt2/{{sample}}/{{sample}}{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
#                 reads = read_ext
#             )
#     output:
#         pic = directory(hic_data_path + "/hic_results/pic/{sample}")
#     params:
#         indir = hic_data_path + "/bowtie_results/bwt2",
#         outdir = hic_data_path
#     log: "logs/hicpro/hicpro_qc_{sample}.log"
#     threads: config['hicpro']['ncpu']
#     shell:
#         """
#         ######################################
#         ## Specific to phoenix for now only ##
#         ######################################
#         ## Load modules
#         module load HiC-Pro/2.9.0-foss-2016b

#         ##Run HiC-pro responding to yes to any interactive requests
#         HiC-Pro \
#           -s quality_checks \
#           -c {input.config} \
#           -i {params.indir} \
#           -o {params.outdir} &> {log}
#         """

# rule hicpro_merge:
#     input:
#         config = hicpro_config,
#         files = hic_data_path + "/hic_results/data/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"
#     output:
#         pairs = hic_data_path + "/hic_results/data/{sample}/{sample}_allValidPairs",
#         stat = hic_data_path + "/hic_results/data/{sample}/{sample}_allValidPairs.mergestat",
#     params:
#         indir = hic_data_path + "/hic_results/data",
#         outdir = hic_data_path
#     log: "logs/hicpro/hicpro_merge_{sample}.log"
#     threads: config['hicpro']['ncpu']
#     shell:
#         """
#         ######################################
#         ## Specific to phoenix for now only ##
#         ######################################
#         ## Load modules
#         module load HiC-Pro/2.9.0-foss-2016b

#         ##Run HiC-pro responding to yes to any interactive requests
#         HiC-Pro \
#           -s merge_persample \
#           -c {input.config} \
#           -i {params.indir} \
#           -o {params.outdir} &> {log}
#         """

# rule build_contact_maps:
#     input:
#         config = hicpro_config,
#         pairs = hic_data_path + "/hic_results/data/{sample}/{sample}_allValidPairs"
#     output:
#         bed = expand([hic_data_path + "/hic_results/matrix/{{sample}}/raw/{bin}/{{sample}}_{bin}_abs.bed"],
#                      bin = bins),
#         mat = expand([hic_data_path + "/hic_results/matrix/{{sample}}/raw/{bin}/{{sample}}_{bin}.matrix"],
#                      bin = bins)
#     params:
#         indir = hic_data_path + "/hic_results/data",
#         outdir = hic_data_path
#     log: "logs/hicpro/build_contact_maps_{sample}.log"
#     threads: config['hicpro']['ncpu']
#     shell:
#         """
#         ######################################
#         ## Specific to phoenix for now only ##
#         ######################################
#         ## Load modules
#         module load HiC-Pro/2.9.0-foss-2016b

#         ##Run HiC-pro responding to yes to any interactive requests
#         HiC-Pro \
#           -s build_contact_maps \
#           -c {input.config} \
#           -i {params.indir} \
#           -o {params.outdir} &> {log}
#         """
