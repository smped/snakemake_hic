rule collect_output_valid_pairs:
    ## This rule collects all final output into the 'output' folder,
    ## enabling safe deletion of all output marked as temp
    input:
        valid_pairs = os.path.join(
            hic_data_path, "hic_results", "data",
            "{sample}/{sample}.allValidPairs"
        )
    output:
        valid_pairs = os.path.join(
            hic_output_path, "allValidPairs", "{sample}.allValidPairs.gz"
        )
    threads: 1
    shell:
        """
        gzip -c {input.valid_pairs} > {output.valid_pairs}
        """

rule collect_output_matrices:
    input:
        matrix = os.path.join(
            hic_data_path, "hic_results", "matrix", "{sample}", "raw", "{bin}",
            "{sample}_{bin}.matrix"
        )
    output:
        matrix = os.path.join(
          hic_output_path, "matrix", "raw", "{bin}", "{sample}_{bin}.matrix.gz"
        )
    threads: 1
    shell:
        """
        gzip -c {input.matrix} > {output.matrix}
        """

rule collect_stats:
    input:
        stat = os.path.join(
            hic_data_path, "hic_results", "stats", "{sample}",
            "{sample}{suffix}"
        )
    output:
        stat = os.path.join(
            hic_output_path, "stats", "{sample}", "{sample}{suffix}"
        )
    threads: 1
    shell:
        """
        cp {input.stat} {output.stat}
        """

rule collect_pics:
    ## Given the problems generating some of these under HiC-Pro, most,
    ## with the exception of the FragmentLength plot can be recreated in a more
    ## informative manner using ggplot directly from the stats output
    ## There is a cheeky copy all line here to cover it though
    input:
        pdf = os.path.join(
            hic_data_path, "hic_results", "pic", "{sample}",
            "plot{file}_{sample}.pdf"
        )
    output:
        pdf = os.path.join(
            hic_output_path, "pic", "{sample}", "plot{file}_{sample}.pdf"
        )
    params:
        in_path = os.path.join(
            hic_data_path, "hic_results", "pic", "{sample}"
        ),
        out_path = os.path.join(hic_output_path, "pic", "{sample}")
    threads: 1
    shell:
        """

        ## The generic copy all
        cp {params.in_path}/*pdf {params.out_path}

        ## And the specific file we know is created & tracked
        cp {input.pdf} {output.pdf}
        """