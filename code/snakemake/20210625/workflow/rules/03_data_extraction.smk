rule get_p_values:
    input:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/associations_raw/{efo_id}.rds")
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/p-values/{efo_id}.txt")
    log:
        os.path.join(config["log_dir"], "get_p_values/{date}/{efo_id}.log")
    container:
        config["R"]
    script:
        "../scripts/get_p_values.R"

rule clump_snps:
    input:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/vcfs/no_dups/{efo_id}.vcf.gz"),
        p_values = os.path.join(config["lts_dir"], "gwasrapidd/{date}/p-values/{efo_id}.txt")
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/plink/clumped/{efo_id}.log")
    log:
        os.path.join(config["log_dir"], "clump_snps/{date}/{efo_id}.log")
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
          --out {params.pref} \
            > {log} 2>&1
        """
# Note: some variants missing from the main dataset, e.g.
# Warning: 'rs199921354' is missing from the main dataset, and is a top variant.
# 1 more top variant ID missing; see log file.

# Note: For some traits, we get the following warning:
# Warning: No significant --clump results.  Skipping.
# In these cases no `*.clumped` file is produced.
# For this reason, we have specified the `*.log` file as the output,
# So that it doesn't cause the snakemake process to fail.

rule get_ref_mafs:
    input:
        all_mafs = os.path.join(config["lts_dir"], "mafs/1kg/20150319/all/all.csv"),
        ref_snps = os.path.join(config["lts_dir"], "gwasrapidd/{date}/plink/clumped/{efo_id}.clumped"),
        # Include log file as input because that was the output of the `clump_snps` rule
        clump_log = os.path.join(config["lts_dir"], "gwasrapidd/{date}/plink/clumped/{efo_id}.log")
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/references/{efo_id}.csv")
    log:
        os.path.join(config["log_dir"], "get_ref_mafs/{date}/{efo_id}.log")
    params:
        percent_interval = config["percent_interval"]
    resources:
        mem_mb = 30000
    container:
        config["R"]
    script:
        "../scripts/get_ref_mafs.R"

rule get_control_snps:
    input:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/references/{efo_id}.csv")
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/control_snps/{efo_id}.txt")
    log:
        os.path.join(config["log_dir"], "get_control_snps/{date}/{efo_id}.log")
    run:
        df = pd.read_csv(input[0])
        df['CONTROL_ID'].to_csv(output[0], header = False, index = False)

rule get_control_vcf_chr:
    input:
        vcf = os.path.join(config["lts_dir"], "vcfs/1kg/20150319/reheaded/{chr}.vcf.gz"),
        snps = os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/control_snps/{efo_id}.txt")
    output:
        os.path.join(config["working_dir"], "vcfs/1kg/20150319/controls/{date}/{efo_id}/by_chr/{chr}.vcf.gz")
    log:
        os.path.join(config["log_dir"], "get_control_vcf_chr/{date}/{efo_id}/{chr}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --include ID=@{input.snps} \
            --output-type z \
            --output-file {output} \
            {input.vcf}
        """

rule merge_control_vcf:
    input:
        expand(os.path.join(config["working_dir"], "vcfs/1kg/20150319/controls/{{date}}/{{efo_id}}/by_chr/{chr}.vcf.gz"),
            chr = CHRS)
    output:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/vcfs/controls/{efo_id}.vcf.gz")
    log:
        os.path.join(config["log_dir"], "merge_control_vcf/{date}/{efo_id}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools concat \
            --output {output.vcf} \
            --output-type z \
            {input}
        """

# Set rule order to avoid `AmbiguousRuleException` as both rules are sending to the `no_dups` directory
ruleorder: merge_control_vcf > get_duplicated_sites

rule get_fst:
    input:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/vcfs/no_dups/{efo_id}.vcf.gz"),
        pop_file = config["local_pop_file"]
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/per_trait/{efo_id}.rds")
    log:
        os.path.join(config["log_dir"], "get_fst/{date}/{efo_id}.log")
    container:
        config["R"]
    script:
        "../scripts/get_fst.R"

rule consolidate_fst:
    input:
        rds = expand(os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/per_trait/{efo_id}.rds"),
            date = DATE_OF_COLLECTION,
            efo_id = EFO_IDS_FILT_PLUS_CONTROL
            )
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/consol/all.rds")
    log:
        os.path.join(config["log_dir"], "consolidate_fst/{date}/all.log")
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/consolidate_fst.R"
