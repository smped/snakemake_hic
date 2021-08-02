rule merge_interaction_matrices:
    input:
        mat = expand(
          hic_data_path + "/hic_results/matrix/{sample}/raw/{{bin}}/{sample}_{{bin}}.matrix",
          sample = samples
        ),
        bed = expand(
          hic_data_path + "/hic_results/matrix/{sample}/raw/{{bin}}/{sample}_{{bin}}_abs.bed",
          sample = samples
        )
    output:
        mat = hic_data_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}.matrix",
        bed = hic_data_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}_abs.bed"
    params:
        samples = samples,
        bin = "{bin}",
        in_path = hic_data_path + "/hic_results/matrix/",
        out_path = hic_data_path + "/hic_results/matrix/merged/raw/"
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