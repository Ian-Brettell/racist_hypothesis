
rule get_fst:
    input:
        vcfs = expand("data/gwasrapidd/{date}/vcfs/original/{efo_id}.vcf.gz",
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

rule clump_snps:
    input:
        vcf = "data/gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz",
        p_values = "data/gwasrapidd/{date}/p-values/{efo_id}.txt"
    output:
        "data/gwasrapidd/{date}/plink/clumped/{efo_id}.clumped"
    params:
        pref = lambda wildcards, output: os.path.splitext(str(output))[0]
    container:
        config["plink1.9"]
    shell:
        """
        plink \
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
