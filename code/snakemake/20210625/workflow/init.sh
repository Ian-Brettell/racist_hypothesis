# For snakemake

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/human_traits_fst
conda activate snakemake_6.7.0
smk_proj="20210625"
snakemake \
  --jobs 5000 \
  --latency-wait 300 \
  --cluster-config code/snakemake/$smk_proj/config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s code/snakemake/$smk_proj/workflow/Snakefile \
  -p

# For R

## Create singularity container on Codon (from `human_traits_fst` repo)
singularity build --remote \
      /hps/nobackup/birney/users/ian/containers/human_traits_fst/R_4.1.0.sif \
      /hps/software/users/birney/ian/repos/human_traits_fst/code/snakemake/20210625/workflow/envs/r_4.1.0/r_4.1.0.def

ssh proxy-codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
# New code:
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  /hps/nobackup/birney/users/ian/containers/human_traits_fst/R_4.1.0.sif

# Old code: was made obsolate when Docker changed its automated build rules
#singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
#                  --bind /hps/software/users/birney/ian/tmp:/tmp \
#                  --bind /hps/software/users/birney/ian/run:/run \
#                  docker://brettellebi/human_traits_fst:R_4.1.0

# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/human_traits_fst/code/snakemake/20210625/workflow/envs/rstudio_server/rsession.conf

ssh -L 8787:hl-codon-44-04:8787 proxy-codon