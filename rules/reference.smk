rule get_reference:
    output: protected(os.path.join(ref_path, ref_fagz))
    params:
        ftp = os.path.join(
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human",
            "release_" + gencode,
            build + "_mapping",
            ref_fagz
        )
    threads: 1
    log: "logs/reference/get_reference.log"
    shell:
        """
        wget \
          -O "{output}" \
          -o {log} \
          {params.ftp}
        """

rule unzip_reference:
    input: rules.get_reference.output
    output: temp(os.path.join(ref_path, ref_fa))
    threads: 1
    shell:
        """
        gunzip -c {input} > {output}
        """

rule bowtie2_index:
    input: rules.unzip_reference.output
    output:
        expand(
            [ref_path + "/bt2/{prefix}.{sub}.bt2"],
            prefix = build + "." + assembly,
            sub = ['1', '2', '3', '4', 'rev.1', 'rev.2']
        )
    conda: "../envs/bowtie2.yml"
    threads: 8
    log: "logs/bowtie2/bowtie2_index.log"
    params:
        prefix = os.path.join(
            ref_path, 
            "bt2", 
            build + "." + assembly
        )
    shell:
        """
        bowtie2-build \
          --threads {threads} \
          -f {input} \
          {params.prefix} &> {log}
        """

rule get_chrom_sizes:
    output: chr_sizes
    params:
        ftp = os.path.join(
            "ftp.ebi.ac.uk/pub/databases/ena/assembly",
            genbank[0:7], 
            genbank[0:10],
            genbank + "_sequence_report.txt"
        ),
    threads: 1
    log: "logs/reference/get_chrom_sizes.log"
    shell:
        """
        # Download the assembly report
        TEMPDIR=$(mktemp -d -t chrXXXXXXXXXX)
        REPORT="assembly_report.txt"
        wget \
            -O "$TEMPDIR/$REPORT" \
            -o {log} \
            {params.ftp}

        # Extract the chrom_sizes
        egrep 'assembled-molecule' "$TEMPDIR/$REPORT" | \
          awk '{{print "chr"$2"\t"$3}}' > {output}

        rm -rf $TEMPDIR
        """