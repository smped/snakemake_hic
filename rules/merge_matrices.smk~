rule merge_interaction_matrices:
    input:
        mat=expand(
          hic_path + "/hic_results/matrix/{sample}/raw/{{bin}}/{sample}_{{bin}}.matrix",
          sample = samples
        ),
        bed=expand(
          hic_path + "/hic_results/matrix/{sample}/raw/{{bin}}/{sample}_{{bin}}_abs.bed",
          sample = samples
        )
    output:
        mat = hic_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}.matrix",
        bed = hic_path + "/hic_results/matrix/merged/raw/{bin}/merged_{bin}_abs.bed"
    conda: "../envs/rtracklayer.yml"
    log: "logs/merge/merge_{bin}.log"
    threads: 4
    shell:
        """
        Rscript --vanilla \
          scripts/merge_matrices.R \
          {input.mat[0]} \
          {input.mat[1]} \
          {input.bed[0]} \
          {input.bed[1]} \
          {output.mat} \
          {output.bed} &> {log}
        """
