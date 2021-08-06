rule initialise_r_project:
    output: rproj
    threads: 1
    shell:
        """
        echo -e "Version: 1.0\n\nRestoreWorkspace: Default\nSaveWorkspace: Default\nAlwaysSaveHistory: Default\n\nEnableCodeIndexing: Yes\nUseSpacesForTab: Yes\nNumSpacesForTab: 2\nEncoding: UTF-8\n\nRnwWeave: knitr\nLaTeX: pdfLaTeX\n\nAutoAppendNewline: Yes\nStripTrailingWhitespace: Yes" > {output}
        """

rule build_site_index:
    input:
        rproj = rproj,
        rulegraph = "rules/rulegraph.dot",
        wflow_yml = "analysis/_site.yml",
        config_yml = "config/config.yml",
        site_yml = "analysis/_site.yml",
        rmd = "analysis/index.Rmd",
        html = expand(
            ["docs/{file}.html"],
            file = ['qc_raw', 'qc_trimmed', 'qc_hic','define_interactions']
        )
    output:
        html = "docs/index.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/build_site_index.log"
    threads: 1
    shell:
        """
        git add analysis/*
        R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
        git add docs/*
        git commit -m 'Updated all'
        """

rule build_raw_qc_report:
    input:
        rproj = rproj,
        fqc = expand(
            ["data/raw/FastQC/{sample}_{reads}_fastqc.zip"],
            reads = ['R1', 'R2'],
            sample = list(df['path'])
        ),
        config_yml = "config/config.yml",
        site_yml = "analysis/_site.yml",        
        rmd = "analysis/qc_raw.Rmd"
    output:
        html = "docs/qc_raw.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/build_raw_qc_report.log"
    threads: 1
    shell:
        """
        R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
        """

rule build_trimmed_qc_report:
    input:
        rproj = rproj,
        fqc = expand(
            ["data/trimmed/FastQC/{sample}_{reads}_fastqc.zip"],
            reads = ['R1', 'R2'],
            sample = list(df['path'])
        ),
        config_yml = "config/config.yml",
        site_yml = "analysis/_site.yml",        
        rmd = "analysis/qc_trimmed.Rmd"
    output:
        html = "docs/qc_trimmed.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/build_trimmed_qc_report.log"
    threads: 1
    shell:
        """
        R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
        """

rule build_hic_qc_report:
    input:
        rproj = rproj,
        frags = fragment_lengths,
        stat =  expand(
            [
                os.path.join(
                    hic_output_path, "stats", "{sample}", "{sample}{suffix}"
                )
            ],
            sample = samples,
            suffix = [
                '.mRSstat', read_ext[0] + ".mmapstat",
                read_ext[1] + ".mmapstat", "_allValidPairs.mergestat"
                ]
            ),
        pic = expand(
                [
                    os.path.join(
                        "docs", "assets", "plotHiCFragmentSize_{sample}.png"
                    )
                ],
                sample = samples
            ),
        config_yml = "config/config.yml",
        site_yml = "analysis/_site.yml",            
        rmd = "analysis/qc_hic.Rmd"
    output:
        html = "docs/qc_hic.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/build_hic_qc_report.log"
    threads: 1
    shell:
        """
        R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
        """

rule define_significant_interactions:
    input:
        rproj = rproj,
        interactions = expand(
            [
                os.path.join(
                    "output", "MaxHiC", "{bin}", "cis_interactions.txt.gz"
                )
            ],
            bin = bins
        ),
        bed = expand(
            [
                os.path.join(
                    hic_output_path, "matrix", "merged_{bin}_abs.bed.gz"
                )
            ],
            bin = bins
        ),
        config_yml = "config/config.yml",
        site_yml = "analysis/_site.yml",            
        rmd = "analysis/define_interactions.Rmd"
    output:
        html = "docs/define_interactions.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/define_interactions.log"
    threads: 8
    shell:
        """
        R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
        """
