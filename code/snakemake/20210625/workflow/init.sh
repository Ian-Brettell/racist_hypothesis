# For snakemake

module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
cd /hps/software/users/birney/ian/repos/human_traits_fst
conda activate snakemake_6.4.1
smk_proj="20210625"
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config code/snakemake/$smk_proj/config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -outdir {cluster.outdir} -o {cluster.outfile} -e {cluster.error}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s code/snakemake/$smk_proj/workflow/Snakefile \
  -p

# For R

bsub -M 20000 -Is """
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  docker://brettellebi/human_traits_fst:R_4.1.0
"""
# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/human_traits_fst/code/snakemake/20210625/workflow/envs/rstudio_server/rsession.conf