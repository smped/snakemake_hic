rule get_reference:
    output:
        fa = temp(os.path.join(ref_path, ref_fa)),
        fagz = os.path.join(ref_path, ref_fagz)
    params:
        genbank = config['ref']['genbank'],
        gencode = config['ref']['gencode'],
        build = config['ref']['build']
    threads: 1
    log: "logs/reference/wget.log"
    shell:
        """
        # Define the URL and download
        URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{params.gencode}/{params.build}_mapping/$(basename {output.fagz})"
        wget \
          -O "{output.fagz}" \
          -o {log} \
          $URL
        gunzip -c "{output.fagz}" > {output.fa}
        """

rule bowtie2_index:
    input: rules.get_reference.output.fa
    output:
        expand([ref_path + "/bt2/{prefix}.{sub}.bt2"],
               prefix = build + "." + assembly,
               sub = ['1', '2', '3', '4', 'rev.1', 'rev.2'] )
    conda: "../envs/bowtie2.yml"
    threads: 8
    log: "logs/bowtie2/bowtie2_index.log"
    params:
        prefix = os.path.join(ref_path, "bt2", build + "." + assembly)
    shell:
        """
        bowtie2-build \
          --threads {threads} \
          -f {input} \
          {params.prefix} &> {log}
        """
