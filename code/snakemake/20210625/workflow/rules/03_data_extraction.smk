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

# 541/3459 traits failed with the following error:
# Error: No variants in VCF file.
# So instead of using the full list of traits, we use the list of traits that both:
# 1. Had >= 1 trait in the VCF
# 2. Had >= significant --clump result (see below)
# ... and therefore produced a .clumped file 
# N passes = 2045, listed in `code/snakemake/20210625/config/20210809_filtered_traits.tsv`
# 541 with no variants + 873 with no significant --clump results + 2045 passes = 3459 total
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

# Note: For 873/3059 traits, we get the following warning:
# Warning: No significant --clump results.  Skipping.
# In these cases no `*.clumped` file is produced.
# For this reason, we have specified the `*.log` file as the output,
# So that it doesn't cause the snakemake process to fail.

rule get_ref_mafs:
    input:
        all_mafs = os.path.join(config["lts_dir"], "mafs/1kg/20201028/all/all.csv"),
        ref_snps = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/plink/clumped/{efo_id}.clumped"),
        # Include log file as input because that was the output of the `clump_snps` rule
        clump_log = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/plink/clumped/{efo_id}.log")
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/references/{efo_id}.csv")
    log:
        os.path.join(config["log_dir"], "get_ref_mafs/{date}/{efo_id}.log")
    params:
        percent_interval = config["percent_interval"], 
        seed = lambda wildcards: filt_traits[filt_traits["EFO_ID"] == wildcards.efo_id]["SEED"]
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
        vcf = os.path.join(config["lts_dir"], "vcfs/1kg/20201028/2504_samples/{chr}.vcf.gz"),
        snps = os.path.join(config["lts_dir"], "gwasrapidd/{date}/controls/control_snps/{efo_id}.txt")
    output:
        os.path.join(config["working_dir"], "vcfs/1kg/20201028/controls/{date}/{efo_id}/by_chr/{chr}.vcf.gz")
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
        expand(os.path.join(config["working_dir"], "vcfs/1kg/20201028/controls/{{date}}/{{efo_id}}/by_chr/{chr}.vcf.gz"),
            chr = CHRS)
    output:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/vcfs/controls/{efo_id}.vcf.gz")
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

rule get_fst_hits:
    input:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/vcfs/no_dups/{efo_id}.vcf.gz"),
        pop_file = config["local_pop_file"]
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/hits/{efo_id}.rds")
    log:
        os.path.join(config["log_dir"], "get_fst/{date}/{efo_id}.log")
    container:
        config["R"]
    script:
        "../scripts/get_fst.R"

rule get_fst_controls:
    input:
        vcf = os.path.join(config["lts_dir"], "gwasrapidd/{date}/high_cov/vcfs/controls/{efo_id}.vcf.gz"),
        pop_file = config["local_pop_file"]
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/controls/{efo_id}.rds")
    log:
        os.path.join(config["log_dir"], "get_fst/{date}/{efo_id}.log")
    container:
        config["R"]
    script:
        "../scripts/get_fst.R"

rule consolidate_fst_hits:
    input:
        rds = expand(os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/hits/{efo_id}.rds"),
            date = DATE_OF_COLLECTION,
            efo_id = EFO_IDS_FILT
            )
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/consol/hits.rds")
    log:
        os.path.join(config["log_dir"], "consolidate_fst/{date}/all.log")
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/consolidate_fst.R"

rule consolidate_fst_controls:
    input:
        rds = expand(os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/controls/{efo_id}.rds"),
            date = DATE_OF_COLLECTION,
            efo_id = EFO_IDS_FILT
            )
    output:
        os.path.join(config["lts_dir"], "gwasrapidd/{date}/pegas/fst/consol/controls.rds")
    log:
        os.path.join(config["log_dir"], "consolidate_fst/{date}/all.log")
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/consolidate_fst.R"
