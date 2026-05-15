cd scattering
snakemake --use-conda --cores all --rerun-incomplete
cd -
cd spectrum
snakemake --use-conda --cores all --rerun-incomplete
cd -