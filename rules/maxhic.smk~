rule get_maxhic:
    output:
        dir = temp(directory("scripts/MaxHiC"))
    threads: 1
    shell:
        """
        wget https://github.com/Rassa-Gvm/MaxHiC/archive/master.zip
        unzip master.zip -d MaxHiC
        mv MaxHiC/MaxHiC-master scripts/MaxHiC
        rm master.zip
        rmdir MaxHiC
        """

rule run_maxhic:
    input:
        maxhic_dir = rules.get_maxhic.output.dir,
        mat = rules.merge_interaction_matrices.output.mat,
        bed = rules.merge_interaction_matrices.output.bed
    output:
        cis = "output/MaxHiC/merged/{bin}/cis_interactions.txt.gz",
        trans = "output/MaxHiC/merged/{bin}/trans_interactions.txt.gz"
    conda: "../envs/maxhic.yml"
    log: "logs/MaxHiC/merged_{bin}_MaxHiC.log"
    threads: 16
    shell:
        """
        ## Given the problems with the raw output from HiC-Pro, we should
        ## delete the *ord.bed* file
        HICDIR=$(dirname {input.mat})
        OUTDIR=$(dirname {output.cis})

        if compgen -G "$HICDIR/*ord.bed" > /dev/null; then
          echo -e "Deleting unnecessary symlink"
          rm $HICDIR/*ord.bed
        fi

        python scripts/MaxHiC/Main.py \
          -t {threads} \
          $HICDIR \
          $OUTDIR &> {log}

        ## Compress the output files
        gzip $OUTDIR/*txt
        """

