rule merge_interaction_matrices:
    input:
        mat = expand(
            os.path.join(
                hic_data_path, "hic_results", "matrix", "{sample}", "raw",
                "{{bin}}", "{sample}_{{bin}}.matrix"
            ),
            sample = samples
        ),
        bed = expand(
            os.path.join(
                hic_data_path, "hic_results", "matrix", "{sample}", "raw",
                "{{bin}}/{sample}_{{bin}}_abs.bed"
            ),
            sample = samples
        )
    output:
        mat = temp(
            os.path.join(
                hic_output_path, "matrix", "raw", "{bin}", "merged_{bin}.matrix"
            )
        ),
        bed = temp(
            os.path.join(
                hic_output_path, "matrix", "raw", "{bin}",
                "merged_{bin}_abs.bed"
            )
        )
    params:
        samples = samples,
        bin = "{bin}",
        in_path = hic_data_path + "/hic_results/matrix/",
        out_path = os.path.join(hic_output_path, "matrix", "raw")
    conda: "../envs/merge_matrices.yml"
    log: "logs/merge/merge_{bin}.log"
    threads: 6
    shell:
        """
        Rscript --vanilla \
            scripts/merge_matrices.R \
            {params.bin} \
            {params.in_path} \
            {params.out_path} \
            {params.samples} &> {log}
        """

rule compress_merged_output:
    input:
        mat = os.path.join(
            hic_output_path, "matrix", "raw", "{bin}", "merged_{bin}.matrix"
        ),
        bed = os.path.join(
            hic_output_path, "matrix", "raw", "{bin}", "merged_{bin}_abs.bed"
        )
    output:
        mat = os.path.join(
            hic_output_path, "matrix", "merged_{bin}.matrix.gz"
        ),
        bed = os.path.join(
            hic_output_path, "matrix", "merged_{bin}_abs.bed.gz"
        )
    threads: 4
    conda: "../envs/pigz.yml"
    shell:
        """
        pigz -p {threads} -c {input.mat} > {output.mat}
        pigz -p {threads} -c {input.bed} > {output.bed}
        """
