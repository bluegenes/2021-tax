"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd

sample_info = pd.read_csv(config['sample_info'], sep = ',')
sample_column = config.get('sample_colname', 'Genome')
SAMPLES = sample_info[sample_column].to_list()
data_dir = config['data_dir']
out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")
thumper_dir = os.path.join(out_dir, "thumper")
basename = config.get('basename', 'tax-basename')
database_dir = 'databases'
search_databases = config.get('search_databases', ['gtdb-rs202-reps'])
alphabet = config.get('alphabet', 'protein')
ksizes = config.get('ksizes', [7,8,9,10,11])
if alphabet in ['nucleotide', 'dna', 'rna']:
    search_type = 'nucleotide'
else:
    search_type = 'protein'

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
        expand(os.path.join(thumper_dir, 'classify', f"{basename}.protein.{alphabet}-k{{ksize}}.classifications.csv"), ksize=ksizes)


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
    input: expand(os.path.join(data_dir, "{accession}.fa.faa"), accession=SAMPLES),
    output: 
        sample_info=os.path.join(out_dir, f"{basename}.protein.thumper.csv")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                sample_name = os.path.basename(inF).rsplit('.fa.faa')[0]
                full_filename = os.path.abspath(str(inF))
                outF.write(f'{sample_name},{full_filename}\n')


rule write_thumper_config:
    input:
        sample_info=os.path.join(out_dir, f"{basename}.{search_type}.thumper.csv")
    output: 
        config= os.path.join(out_dir, f"{basename}.{alphabet}.thumper.conf") 
    params:
        ksize=ksizes,
    run:
        with open(str(output), 'w') as out:
            out.write(f"strict: 1\n")
            out.write(f"force: 0\n")
            out.write(f"input_type: protein\n")
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
         expand(os.path.join(thumper_dir, 'classify', f"{basename}.protein.{alphabet}-k{{ksize}}.classifications.csv"),ksize=ksizes)
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "thumper_classify", f"{basename}.protein.{alphabet}.log")
    benchmark: os.path.join(logs_dir, "thumper_classify", f"{basename}.protein.{alphabet}.benchmark")
    conda: "conf/env/thumper.yml"
    shell:
        """
        thumper run {input.config} --jobs {threads} genome_classify --nolock
        """
