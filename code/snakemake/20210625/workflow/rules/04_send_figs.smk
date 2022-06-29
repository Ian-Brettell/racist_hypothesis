# Send plots to be used for the paper to the Google Drive
rule send_plots_to_google_drive:
    input:
        final = expand(os.path.join(
            config["repo_dir"],
            "docs/plots/{clump_params}/{clump_params}_20220314_final.{ext}"
            ),
                clump_params = config["clump_r2_kb"],
                ext = ["png", "pdf"]
        ),
        snp_counts = expand(os.path.join(
            config["repo_dir"],
            "docs/plots/0.1_1000/0.1_1000_20220314_snp_counts.{ext}"
            ),
                ext = ["png", "pdf"]
        ),
        ecdf_all = expand(os.path.join(
            config["repo_dir"],
            "docs/plots/0.1_1000/0.1_1000_20220314_ecdf_all_faceted_with_controls_d_rank_long.{ext}"
            ),
                ext = ["png", "pdf"]
        ),
    output:
        touch(
            os.path.join(
                config["working_dir"],
                "logs/send_plots_to_google_drive/all.done"
            )
        )
    log:
        os.path.join(
            config["working_dir"],
            "logs/send_plots_to_google_drive/all.log"
        ),
    params:
        drive_dir = config["google_drive_dir"]
    conda:
        config["rclone_env"]
    resources:
        mem_mb = 1000
    shell:
        """
        for i in $(echo {input}); do \
            rclone copy $i {params.drive_dir}/ ; \
        done
        """
