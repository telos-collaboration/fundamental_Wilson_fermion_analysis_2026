cd scattering
snakemake --use-conda --cores 8 --forceall --rerun-incomplete
cd -
cd spectrum
snakemake --use-conda --cores 8 --forceall --rerun-incomplete
cd -