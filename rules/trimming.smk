rule adapter_removal:
    input:
        r1 = "data/raw/fastq/{sample}/{file}" + read_ext[0] + suffix,
        r2 = "data/raw/fastq/{sample}/{file}" + read_ext[1] + suffix
    output:
        t1 = temp(trim_path + "/{sample}/{file}" + read_ext[0] + suffix),
        t2 = temp(trim_path + "/{sample}/{file}" + read_ext[1] + suffix),
        log = "data/trimmed/logs/{sample}/{file}.settings"
    conda:
        "../envs/adapterremoval.yml"
    params:
        adapter1 = config['trimming']['adapter1'],
        adapter2 = config['trimming']['adapter2'],
        minlength = config['trimming']['minlength'],
        minqual = config['trimming']['minqual'],
        maxns = config['trimming']['maxns']
    threads: 4
    log:
        "logs/adapterremoval/{sample}/{file}.log"
    shell:
        """
        AdapterRemoval \
            --adapter1 {params.adapter1} \
            --adapter2 {params.adapter2} \
            --file1 {input.r1} \
            --file2 {input.r2} \
            --threads {threads} \
            --gzip \
            --maxns {params.maxns} \
            --trimqualities \
            --minquality {params.minqual} \
            --minlength {params.minlength} \
            --output1 {output.t1} \
            --output2 {output.t2} \
            --discarded /dev/null \
            --singleton /dev/null \
            --settings {output.log} &> {log}
        """
