rule install_maxhic:
    output:
        exe = temp("scripts/MaxHiC/Main.py"),
        md = temp("scripts/MaxHiC/README.md"),
        general = temp(directory("scripts/MaxHiC/General")),
        capture = temp(directory("scripts/MaxHiC/Capture"))
    log: "logs/MaxHiC/install.log"
    threads: 1
    shell:
        """
        git clone https://github.com/Rassa-Gvm/MaxHiC.git scripts/MaxHiC
        rm -rf scripts/MaxHiC/Sample_Inputs
        """

rule run_maxhic:
    input:
        maxhic_exe = ancient(rules.install_maxhic.output.exe),
        mat = os.path.join(
            hic_output_path, "matrix", "raw", "{bin}", "merged_{bin}.matrix"
        ),
        bed = os.path.join(
            hic_output_path, "matrix", "raw", "{bin}", "merged_{bin}_abs.bed"
        )
    output:
        cis = "output/MaxHiC/{bin}/cis_interactions.txt.gz",
        trans = "output/MaxHiC/{bin}/trans_interactions.txt.gz",
        model_params = directory("output/MaxHiC/{bin}/ModelParameters")
    params:
        input_path = os.path.join(hic_output_path, "matrix", "raw", "{bin}"),
        output_path = os.path.join("output", "MaxHiC", "{bin}")
    conda: "../envs/maxhic.yml"
    log: "logs/MaxHiC/merged_{bin}_MaxHiC.log"
    threads: 16
    shell:
        """
        ## Given the problems with the raw output from HiC-Pro, we should
        ## delete any *ord.bed* files that exist. They seem to have been 
        ## excluded from HiC-Pro v3
        if compgen -G "{params.input_path}/*ord.bed" > /dev/null; then
          echo -e "Deleting unnecessary symlink"
          rm {params.input_path}/*ord.bed
        fi

        python {input.maxhic_exe} \
          -t {threads} \
          {params.input_path} \
          {params.output_path} &> {log}

        ## Compress the output files
        pigz -p {threads} {params.output_path}/*txt
        """
