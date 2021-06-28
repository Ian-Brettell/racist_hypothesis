
rule get_fst:
    input:
        vcfs = expand("data/gwasrapidd/{date}/vcfs/{efo_id}.vcf.gz",
                        date = DATE_OF_COLLECTION,
                        efo_id = EFO_IDS
                    ),
        pop_file = config["local_pop_file"]
    output:
        "data/gwasrapidd/{date}/pegas/fst/fst.rds"
    container:
        config["R"]
    script:
        "../scripts/get_fst.R"

rule get_p_values:
    input:
        "data/gwasrapidd/{date}/associations_raw/{efo_id}.rds"
    output:
        "data/gwasrapidd/{date}/p-values/{efo_id}.txt"
    container:
        config["R"]
    script:
        "../scripts/get_p_values.R"

rule clump_snps:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/{efo_id}.vcf.gz",
        p_values = "data/gwasrapidd/{date}/p-values/{efo_id}.txt"
    output:
        "data/gwasrapidd/{date}/plink/clumped/{efo_id}.log"
    params:
        pref = lambda wildcards, output: os.path.splitext(str(output))[0]
    container:
        config["plink1.9"]
    shell:
        """
        /usr/bin/plink1.9 \
          --vcf {input.vcf} \
          --clump {input.p_values} \
          --clump-p1 0.00000001 \
          --clump-p2 0.00000001 \
          --clump-r2 {config[clump_r2]} \
          --clump-kb {config[clump_kb]} \
          --out {params.pref}     
        """
# Note: some variants missing from the main dataset, e.g.
# Warning: 'rs199921354' is missing from the main dataset, and is a top variant.
# 1 more top variant ID missing; see log file.

# Note: For some traits, we get the following warning:
# Warning: No significant --clump results.  Skipping.
# In these cases no `*.clumped` file is produced.
# For this reason, we have specified the `*.log` file as the output,
# So that it doesn't cause the snakemake process to fail.
