"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd

sample_info = pd.read_csv("raw/tmp-1000r-mgnify-genomes-all_metadata.tsv", sep = '\t')
SAMPLES = sample_info['Genome'].to_list()
data_dir = "raw/v1.0"
out_dir = 'outputs'
logs_dir = os.path.join(out_dir, "logs")
thumper_dir = os.path.join(out_dir, "thumper")
basename = '1000r-mgnify-genomes'
database_dir = 'databases'
search_databases = []
alphabet = 'protein'
ksizes = [7,8,9,10,11]

# steps
# 1. read sample csv
# 2. prodigal fna --> faa
# 3. thumper classify (protein input)
# 4. compare!
#  - compare vs "true" in original csv
#  - compare vs gtdb-tk


rule all:
    input: os.path.join(out_dir, f"{basename}.{alphabet}.thumper.conf")


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


rule write_sample_info:
    input: expand(os.path.join(data_dir, "proteins_v1.0", "{accession}.prodigal.faa"), accession=SAMPLES),
    output: 
        sample_info=os.path.join(out_dir, "{basename}.prodigal.thumper.csv")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                sample_name = os.path.basename(inF).rsplit('.prodigal.faa')[0]
                full_filename = os.path.abspath(str(inF))
                outF.write(f'{sample_name},{full_filename}\n')


rule write_thumper_config:
    input:
        sample_info = rules.write_sample_info.output.sample_info,
    output: 
        config= os.path.join(out_dir, "{basename}.{alphabet}.thumper.conf") 
    params:
        metagenome_trim_memory=config.get("metagenome_trim_memory", "1e9"),
        ksize=ksizes
    run:
        with open(str(output), 'w') as out:
            out.write(f"strict: 1\n")
            out.write(f"force: 0\n")
            out.write(f"input_type: protein\n")
            out.write(f"outdir: {thumper_dir}\n")
            out.write(f"basename: {basename}\n")
            out.write(f"sample_info: {input.sample_info}\n")
            out.write(f"database_dir: {database_dir}\n")
            out.write(f"search_databases:\n")
            for sd in search_databases:
                out.write(f"  - {sd}\n")
            out.write(f"alphabet: {wildcards.alphabet}\n")
            out.write(f"ksize:\n")
            for k in params.ksize:
                out.write(f"  - {k}\n")


rule thumper_classify:
    input: rules.write_thumper_config.output.config
    output:
         os.path.join(thumper_dir, 'classify', f"{basename}.protein.{alphabet}-k{{ksize}}.classifications.csv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "thumper_classify", f"{basename}.protein.{alphabet}-k{{ksize}}.log")
    benchmark: os.path.join(logs_dir, "thumper_classify", f"{basename}.protein.{alphabet}-k{{ksize}}.benchmark")
    conda: "conf/env/thumper.yml"
    shell:
        """
        thumper genome_classify {input} -j {threads} 
        """
