rule find_restriction_fragments:
    input: ancient(os.path.join(ref_path, ref_fa))
    output:
        rs = temp(rs_frags)
    params:
        restriction_site = config['hicpro']['restriction_site'],
        script = os.path.join(
          os.path.dirname(hic_path), 
          "utils", 
          "digest_genome.py"
        )
    threads: 1
    conda: "../envs/hicpro.yml"
    shell:
        """
        # Run the python script
        python {params.script} \
          -r {params.restriction_site} \
          -o {output.rs} \
          {input}
        """

rule get_fragment_lengths:
    input: rs_frags
    output: fragment_lengths
    threads: 1
    shell:
        """
        awk '{{print($3 - $2)}}' {input} | \
          sort -n | \
          uniq -c > {output}
        """

rule zip_fragment_bedfile:
    input: rules.find_restriction_fragments.output.rs
    output: rs_frags + ".gz"
    threads: 4
    conda: "../envs/pigz.yml"
    shell:
        """
        pigz -p {threads} -c {input} > {output}
        """

rule make_hicpro_config:
    input:
        idx = ancient(
          expand(
            [ref_path + "/bt2/{prefix}.{sub}.bt2"],
            prefix = build + "." + assembly,
            sub = ['1', '2', '3', '4', 'rev.1', 'rev.2']
          )
        ),
        rs = rs_frags,
        chr_sizes = chr_sizes,
	config = "config/config.yml"
    output:
        hicpro_config
    params:
        template = os.path.join(
          os.path.dirname(
            os.path.dirname(hic_path)
          ), 
          "config-hicpro.txt"
        ),
        idx_root = os.path.join(
          ref_path, "bt2"
        )
    conda: "../envs/stringr.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla \
          scripts/write_hicpro_config.R \
          {params.idx_root} \
          {input.chr_sizes} \
          {input.rs} \
          {params.template} \
          {output}
        """

rule hicpro_mapping:
    input:
        config = hicpro_config,
        files = expand(
          [trim_path + "/{sample_path}{read_ext}{suffix}"],
          sample_path = df['path'],
          read_ext = read_ext,
          suffix = suffix
        ),
        idx = ancient(
          expand(
            [ref_path + "/bt2/{prefix}.{sub}.bt2"],
            prefix = build + "." + assembly,
            sub = ['1', '2', '3', '4', 'rev.1', 'rev.2']
          )
        )
    output:
        bwt2 = temp(
            expand(
                ["{path}/{sample_path}{reads}_" + build + "." + assembly + "{suffix}"],
                path = os.path.join(bowtie_data_path, "bwt2"),
                sample_path = df['path'],
                reads = read_ext,
                suffix = ['.bwt2merged.bam', '.mapstat']
            )
        ),
        bwt2_global = temp(
          expand(
            ["{path}/{sample_path}{reads}_{meta}.{suffix}"],
            path = os.path.join(bowtie_data_path, "bwt2_global"),
            sample_path = df['path'],
            reads = read_ext,
            meta = build + "." + assembly + ".bwt2glob",
            suffix = ['bam', 'unmap.fastq']
          )
        ),
        bwt2_local = temp(
          expand(
            ["{path}/{sample_path}{reads}_{meta}_{suffix}"],
            path = os.path.join(bowtie_data_path, "bwt2_local"),
            sample_path = df['path'],
            reads = read_ext,
            meta = build + "." + assembly + ".bwt2glob.unmap",
            suffix = ['bwt2loc.bam', 'trimmed.fastq']
          )
        )
    params:
        outdir = hic_data_path,
        indir = trim_path
    threads: config['hicpro']['ncpu']
    conda: "../envs/hicpro.yml"
    shell:
        """
        ## Remove any existing data as leaving this here causes HicPro to
        ## make an interactive request. Piping `yes` into HicPro may be the
        ## source of some current problems
        if [[ -d {params.outdir} ]]; then
          rm -rf {params.outdir}
        fi

        ## Run HiC-pro
        HiC-Pro \
          -s mapping \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} 
        """

rule hicpro_qc:
    input:
        config = hicpro_config,
        bwt2 = expand(
            ["{path}/{sample_path}{reads}_" + build + "." + assembly + "{suffix}"],
            path = os.path.join(bowtie_data_path, "bwt2"),
            sample_path = df['path'],
            reads = read_ext,
            suffix = ['.bwt2merged.bam', '.mapstat']
        ),
        pairs = expand(
            [
                os.path.join(
                    hic_data_path, "hic_results", "data", "{sample}", 
                    "{sample}.allValidPairs"
                )
            ],
            sample = samples
        ),
	      stat = expand(
	          [
            	os.path.join(
                    hic_data_path, "hic_results", "stats", "{sample}",
                    "{sample}{suffix}"
                )
            ],
            sample = samples,
            suffix = ['.mRSstat', read_ext[0] + ".mmapstat", read_ext[1] + ".mmapstat", "_allValidPairs.mergestat", ".mpairstat"]
        )
    output:
        temp(
            expand(
                [
                    os.path.join(
                        hic_data_path, "hic_results", "pic", "{sample}",
                        "plot{file}_{sample}.pdf"
                    )
                ],
                sample = samples,
                file = ['HiCFragment', 'MappingPairing', 'Mapping', 'HiCContactRanges', 'HiCFragmentSize']
            )
        )
    params:
        indir = os.path.join(bowtie_data_path, "bwt2"),
        outdir = hic_data_path
    threads: 1
    conda: "../envs/hicpro.yml"
    shell:
        """
        HiC-Pro \
          -s quality_checks \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} 
        """
        
rule hicpro_proc:
    input:
        config = hicpro_config,
        bwt2 = expand(
            ["{path}/{sample_path}{reads}_" + build + "." + assembly + "{suffix}"],
            path = os.path.join(bowtie_data_path, "bwt2"),
            sample_path = df['path'],
            reads = read_ext,
            suffix = ['.bwt2merged.bam', '.mapstat']
        )
    output:
        bam = temp(
            expand(
                [bowtie_data_path + "/bwt2/{sample_path}_{meta}.{suffix}"],
                meta = build + "." + assembly + ".bwt2pairs",
                suffix = ['bam', 'pairstat'],
                sample_path = df['path']
              )
        ),
        pairs = temp(
            expand(
                [hic_data_path + "/hic_results/data/{sample_path}_{meta}.{suffix}"],
                meta = build + "." + assembly + ".bwt2pairs",
                suffix = ['DEPairs', 'DumpPairs', 'FiltPairs', 'REPairs', 'RSstat', 'SCPairs', 'SinglePairs', 'validPairs'],
                sample_path = df['path']
            )
        )
    params:
        indir = os.path.join(bowtie_data_path, "bwt2"),
        outdir = hic_data_path
    threads: config['hicpro']['ncpu']
    conda: "../envs/hicpro.yml"
    shell:
        """
        HiC-Pro \
          -s proc_hic \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} 
        """

rule hicpro_merge:
    input:
        config = hicpro_config,
        files = expand(
            [hic_data_path + "/hic_results/data/{sample_path}_{meta}.{suffix}"],
            meta = build + "." + assembly + ".bwt2pairs",
            suffix = ['DEPairs', 'DumpPairs', 'FiltPairs', 'REPairs', 'RSstat', 'SCPairs', 'SinglePairs', 'validPairs'],
            sample_path = df['path']
            )
    output:
        pairs = temp(
            expand(
                [
                    os.path.join(
                        hic_data_path, "hic_results", "data", "{sample}", 
                        "{sample}.allValidPairs"
                    )
                ],
                sample = samples
            )
        ),
        stat = temp(
            expand(
                [
                    os.path.join(
                        hic_data_path, "hic_results", "stats", "{sample}",
                        "{sample}{suffix}"
                    )
                ],
                sample = samples,
                suffix = ['.mRSstat', read_ext[0] + ".mmapstat", read_ext[1] + ".mmapstat", "_allValidPairs.mergestat", ".mpairstat"]
            )
        )
    params:
        indir = hic_data_path + "/hic_results/data",
        outdir = hic_data_path
    threads: config['hicpro']['ncpu']
    conda: "../envs/hicpro.yml"
    shell:
        """
        HiC-Pro \
          -s merge_persample \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} 
        """

rule build_contact_maps:
    input:
        config = hicpro_config,
        pairs = expand(
          [hic_data_path + "/hic_results/data/{sample}/{sample}.allValidPairs"],
          sample = samples
        )
    output:
      matrix = temp(
          expand(
            [
              os.path.join(
                hic_data_path, "hic_results", "matrix", "{sample}", "raw",
                "{bin}", "{sample}_{bin}{suffix}"
              )
            ],
          sample = samples,
          bin = bins,
          suffix = ['.matrix', '_abs.bed']
        )
      )
    params:
        indir = hic_data_path + "/hic_results/data",
        outdir = hic_data_path
    threads: config['hicpro']['ncpu']
    conda: "../envs/hicpro.yml"
    shell:
        """
        HiC-Pro \
          -s build_contact_maps \
          -c {input.config} \
          -i {params.indir} \
          -o {params.outdir} 
        """
