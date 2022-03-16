"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd

configfile: "conf/mgnify-1000r.protein.gtdb-rs202.conf"

sample_info = pd.read_csv(config['sample_info'], sep = '\t')
sample_column = config.get('sample_colname', 'Genome')
lineage_column = config.get('sample_colname', 'Lineage')

SAMPLES = sample_info[sample_column].to_list()
data_dir = config['data_dir']
out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")
thumper_dirname= config.get("thumper_dir", "thumper")
reports_dirname= config.get("reports_dir", "reports")
thumper_dir = os.path.join(out_dir, thumper_dirname)
reports_dir = os.path.join(out_dir, reports_dirname)
basename = config.get('basename', 'tax-basename')
database_dir = 'databases'

search_databases = config.get('search_databases', ['gtdb-rs202'])
input_type = config.get('input_type', 'protein')
alphabet = config.get('alphabet', 'protein')
ksizes = config.get('ksizes', [7,8,9,10,11])
sketch_type= "nucleotide" if alphabet in ['nucleotide','dna', 'rna' ] else "protein"
threshold_bp = config.get('threshold_bp', '50000')

additional_db_info = config.get('thumper_user_database_info')
# steps
# 1. read sample csv
# 2. prodigal fna --> faa
# 3. thumper classify (protein input)
# 4. compare!
#  - compare vs "true" in original csv
#  - compare vs gtdb-tk


rule all:
    input: 
        os.path.join(out_dir, f"{basename}.{alphabet}.thumper.conf"),
        expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"), ksize=ksizes),
        os.path.join(out_dir, "gtdbtk-classify", f"{basename}.nucl-gtdbtk.bac120.summary.tsv"),
        os.path.join(out_dir, "gtdbtk-classify", f"{basename}.nucl-gtdbtk.ar122.summary.tsv"),
        os.path.join(out_dir, "gtdbtk-classify", f"{basename}.gatktk_lineages.csv"),
        expand(os.path.join(reports_dir, f'{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classification-report.csv'),ksize=ksizes),
        os.path.join(reports_dir, f'{basename}.{sketch_type}.classification-report.csv'),
        os.path.join(reports_dir, f'{basename}.{sketch_type}.classification-summaries.csv'),

rule prodigal_genomes:
    input: os.path.join(data_dir, "{accession}.fa")
    output:
        proteins=os.path.join(data_dir, "proteins_v1.0", "{accession}.prodigal.faa"),
        genes=os.path.join(data_dir, "proteins_v1.0","{accession}.prodigal.fna")
    conda: "conf/env/prodigal-env.yml"
    log: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.log")
    benchmark: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.benchmark")
    shell:
        """
         prodigal -i {input} -a {output.proteins} -o {output.genes} 2> {log}
        """


rule write_nucl_sample_info:
    input: expand(os.path.join(data_dir, "{accession}.fa"), accession=SAMPLES),
    output:
        sample_info=os.path.join(out_dir, f"{basename}.nucleotide.thumper.csv")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                sample_name = os.path.basename(inF).rsplit('.fa')[0]
                full_filename = os.path.abspath(str(inF))
                outF.write(f'{sample_name},{full_filename}\n')


rule write_protein_sample_info:
    input: expand(os.path.join(data_dir, "proteins_v1.0", "{accession}.prodigal.faa"), accession=SAMPLES),
    output: 
        sample_info=os.path.join(out_dir, f"{basename}.protein.thumper.csv")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                sample_name = os.path.basename(inF).rsplit('.prodigal.faa')[0]
                full_filename = os.path.abspath(str(inF))
                outF.write(f'{sample_name},{full_filename}\n')


rule write_thumper_config:
    input:
        sample_info=os.path.join(out_dir, f"{basename}.{alphabet}.thumper.csv")
    output: 
        config= os.path.join(out_dir, f"{basename}.{alphabet}.thumper.conf") 
    params:
        ksize=ksizes,
    run:
        with open(str(output), 'w') as out:
            out.write(f"strict: 1\n")
            out.write(f"force: 0\n")
            out.write(f"input_type: {input_type}\n")
            out.write(f"data_dir: {data_dir}\n")
            out.write(f"output_dir: {thumper_dir}\n")
            out.write(f"basename: {basename}\n")
            out.write(f"sample_info: {input.sample_info}\n")
            out.write(f"database_dir: {database_dir}\n")
            out.write(f"sourmash_database_threshold_bp: {threshold_bp}\n")
            out.write(f"search_databases:\n")
            for sd in search_databases:
                out.write(f"  - {sd}\n")
            if additional_db_info is not None:
                out.write(f"user_database_info: {additional_db_info}\n")
            out.write(f"alphabet: {alphabet}\n")
            out.write(f"ksize:\n")
            for k in params.ksize:
                out.write(f"  - {k}\n")


rule thumper_classify:
    input: 
        config=rules.write_thumper_config.output.config
    output:
         expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"),ksize=ksizes)
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "thumper_classify", f"{basename}.{sketch_type}.{alphabet}.log")
    benchmark: os.path.join(logs_dir, "thumper_classify", f"{basename}.{sketch_type}.{alphabet}.benchmark")
    conda: "conf/env/thumper.yml"
    shell:
        """
        thumper run {input.config} --jobs {threads} genome_classify --nolock
        """


rule write_gtdbtk_batchfile:
    input: expand(os.path.join(data_dir, "{accession}.fa"), accession=SAMPLES),
    output:
        sample_info=os.path.join(out_dir, f"{basename}.nucleotide.batchfile.tsv")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                sample_name = os.path.basename(inF).rsplit('.fa')[0]
                full_filename = os.path.abspath(str(inF))
                outF.write(f'{full_filename}\t{sample_name}\n')


rule gtdbtk_classify_from_nucleotide:
    input:
        sample_info=os.path.join(out_dir, "{basename}.nucleotide.batchfile.tsv")
    output:
        os.path.join(out_dir, "gtdbtk-classify", "{basename}.nucl-gtdbtk.bac120.summary.tsv"),
        os.path.join(out_dir, "gtdbtk-classify", "{basename}.nucl-gtdbtk.ar122.summary.tsv"),
    params:
        outD= os.path.join(out_dir, "gtdbtk-classify"),
        scratch=config.get('scratch_dir', '/scratch/ntpierce'),
    threads: 64
    resources:
        mem_mb=lambda wildcards, attempt: attempt *210000,
        runtime=1200,
    log: os.path.join(logs_dir, "gtdbtk", "{basename}.gtdbtk.log")
    benchmark: os.path.join(logs_dir, "gtdbtk", "{basename}.gtdbtk.benchmark")
    conda: "conf/env/gtdbtk.yml"
    shell:
        """
        export GTDBTK_DATA_PATH=/group/ctbrowngrp/gtdb/gtdbtk/release202
        gtdbtk classify_wf --batchfile {input} --out_dir {params.outD} -x ".fa" \
                           --cpus {threads} --scratch_dir {params.scratch} \
                           --prefix "{wildcards.basename}.nucl-gtdbtk" --force \
                           2> {log}
        """


rule gtdbtk_to_reference_lineage_csv:
    input: 
        bac120=os.path.join(out_dir, "gtdbtk-classify", "{basename}.nucl-gtdbtk.bac120.summary.tsv"),
        ar122=os.path.join(out_dir, "gtdbtk-classify", "{basename}.nucl-gtdbtk.ar122.summary.tsv"),
    output:
        os.path.join(out_dir, "gtdbtk-classify", "{basename}.gatktk_lineages.csv")
    run:
        bac120 = pd.read_csv(str(input.bac120), sep = '\t')
        ar122 = pd.read_csv(str(input.ar122), sep = '\t')
        bac120 = bac120.append(ar122, ignore_index=True)
        bac120.rename(columns={'user_genome': 'genome', 'classification': 'lineage'}, inplace=True)
        bac120[['genome', 'lineage']].to_csv(str(output), index=False, sep=',')


rule write_reference_lineage_csv:
    input: config['sample_info']
    output:
        os.path.join(out_dir, "{basename}.reference_lineages.csv")
    run:
        ref_lineages = sample_info[[sample_column, lineage_column]]
        ref_lineages.rename(columns={sample_column: 'genome', lineage_column:'lineage'}, inplace=True)
        ref_lineages.to_csv(str(output), index=False)


rule assess_classification_each_ksize:
    input:
        ref_lin = os.path.join(out_dir, "gtdbtk-classify", f"{basename}.gatktk_lineages.csv"),
        tax_csvs = os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"),
    output:
        classif_details = os.path.join(reports_dir, f'{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classification-report.csv'),
        classify_fsummaries = os.path.join(reports_dir, f'{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classification-summaries.csv')
    log: os.path.join(logs_dir, "assess-classification", reports_dirname, f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.log")
    benchmark: os.path.join(logs_dir, "assess-classification", reports_dirname, f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.benchmark")
    #conda: "conf/env/thumper.yml"
    conda: "sourmash-dev.yml"
    shell:
        """
        python assess-classification.py --tax-genome {input.tax_csvs} \
                                        --ref-annotations {input.ref_lin} \
                                         --output-csv {output.classif_details} \
                                         --report-csv {output.classify_fsummaries} \
                                         --reference-type 'gtdbtk' \
                                         2> {log}
        """


rule assess_classification_all:
    input:
        ref_lin = os.path.join(out_dir, "gtdbtk-classify", f"{basename}.gatktk_lineages.csv"),
        tax_csvs = expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"), ksize=ksizes)
    output:
        classif_details = os.path.join(reports_dir, f'{basename}.{sketch_type}.classification-report.csv'),
        classify_fsummaries = os.path.join(reports_dir, f'{basename}.{sketch_type}.classification-summaries.csv')
    log: os.path.join(logs_dir, "assess-classification", reports_dirname, f"{basename}.{sketch_type}.{alphabet}.log")
    benchmark: os.path.join(logs_dir, "assess-classification", reports_dirname, f"{basename}.{sketch_type}.{alphabet}.benchmark")
    #conda: "conf/env/thumper.yml"
    conda: "sourmash-dev.yml"
    shell:
        """
        python assess-classification.py --tax-genome {input.tax_csvs} \
                                        --ref-annotations {input.ref_lin} \
                                         --output-csv {output.classif_details} \
                                         --report-csv {output.classify_fsummaries} \
                                         --reference-type 'gtdbtk' \
                                         2> {log}
        """
