"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd

sample_info = pd.read_csv(config['sample_info'], sep = '\t')
sample_column = config.get('sample_colname', 'Genome')
lineage_column = config.get('sample_colname', 'Lineage')

SAMPLES = sample_info[sample_column].to_list()
data_dir = config['data_dir']
out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")
thumper_dir = os.path.join(out_dir, "thumper")
basename = config.get('basename', 'tax-basename')
database_dir = 'databases'
search_databases = config.get('search_databases', ['gtdb-rs202-reps'])
input_type = config.get('input_type', 'protein')
alphabet = config.get('alphabet', 'protein')
ksizes = config.get('ksizes', [7,8,9,10,11])
sketch_type= "nucleotide" if alphabet in ['nucleotide','dna', 'rna' ] else "protein"
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
        classif_details = os.path.join(out_dir, 'reports', f'{basename}.{sketch_type}.classification-report.csv'),
        classify_fsummaries = os.path.join(out_dir, 'reports', f'{basename}.{sketch_type}.classification-summaries.csv')


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
            out.write(f"search_databases:\n")
            for sd in search_databases:
                out.write(f"  - {sd}\n")
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

rule write_reference_lineage_csv:
    input: config['sample_info']
    output:
        os.path.join(out_dir, f"{basename}.reference_lineages.csv")
    run:
        ref_lineages = sample_info[[sample_column, lineage_column]]
        ref_lineages.rename(columns={sample_column: 'genome', lineage_column:'lineage'}, inplace=True)
        ref_lineages.to_csv(str(output), index=False)


rule assess_classification:
    input:
        ref_lin = rules.write_reference_lineage_csv.output,
        tax_csvs = expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"), ksize=ksizes)
    output:
        classif_details = os.path.join(out_dir, 'reports', f'{basename}.{sketch_type}.classification-report.csv'),
        classify_fsummaries = os.path.join(out_dir, 'reports', f'{basename}.{sketch_type}.classification-summaries.csv')
    log: os.path.join(logs_dir, "assess-classification", f"{basename}.{sketch_type}.{alphabet}.log")
    benchmark: os.path.join(logs_dir, "assess-classification", f"{basename}.{sketch_type}.{alphabet}.benchmark")
    #conda: "conf/env/thumper.yml"
    conda: "sourmash-dev.yml"
    shell:
        """
        python assess-classification.py --tax-genome {input.tax_csvs} \
                                        --ref-annotations {input.ref_lin} \
                                         --output-csv {output.classif_details} \
                                         --report-csv {output.classify_fsummaries} \
                                         2> {log}
        """
