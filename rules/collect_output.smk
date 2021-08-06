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
    ## Given the problems generating some of these under HiC-Pro, most,
    ## with the exception of the FragmentLength plot can be recreated in a more
    ## informative manner using ggplot directly from the stats output
    ## There is a cheeky copy all line here to cover it though
    input:
        stat = os.path.join(
            hic_data_path, "hic_results", "stats", "{sample}",
            "{sample}{suffix}"
        )
    output:
        stat = os.path.join(
            hic_output_path, "stats", "{sample}", "{sample}{suffix}"
        )
    params:
        in_path = os.path.join(
            hic_data_path, "hic_results", "stats", "{sample}"
        ),
        out_path = os.path.join(hic_output_path, "stats", "{sample}")
    threads: 1
    shell:
        """
        ## The generic copy all. Clearly this will repeat every time though
        # cp {params.in_path}/* {params.out_path}

        ## Copy the specific files
        cp {input.stat} {output.stat}
        """

rule convert_fragment_size_plot:
    input: 
        pdf = os.path.join(
            hic_data_path, "hic_results", "pic", "{sample}",
            "plotHiCFragmentSize_{sample}.pdf"
        )
    output:
        png = os.path.join(
            "docs", "assets", "plotHiCFragmentSize_{sample}.png"
        )
    params:
        density = 150,
        quality = 90
    threads: 1
    conda: "../envs/imagemagick.yml"
    shell:
        """
        convert \
            -density {params.density} \
            {input.pdf}[0] \
            -quality {params.quality} \
            {output.png}
        """
