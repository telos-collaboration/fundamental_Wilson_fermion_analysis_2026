# Resonant scattering in two-flavored Sp(4) lattice gauge theories&mdash;Analysis workflow

Code: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20638262.svg)](https://doi.org/10.5281/zenodo.20638262)
Data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20111459.svg)](https://doi.org/10.5281/zenodo.20111459)

The workflow in this repository performs
the analyses presented in the paper
[Resonant scattering in two-flavored Sp(4) lattice gauge theories&mdash;Analysis workflow][paper].

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda


## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/telos-collaboration/fundamental_Wilson_fermion_analysis_2026
   cd fundamental_Wilson_fermion_analysis_2026
   ```

3. Download the archive containing the raw log files named `raw_data.tar` and untar it in this directory.

4. The workflow is run using Snakemake:

    ```shellsession
    snakemake --cores 1 --use-conda
    ```

where the number `1` may be replaced by the number of CPU cores you wish to allocate to the computation.

Snakemake will automatically download and install all required Python packages. This requires an Internet connection; if you are running in an HPC environment where you would need to run the workflow without Internet access, details on how to preinstall the environment can be found in the Snakemake documentation.


## Output

Output plots, tables, and definitions are placed into the `assets` directory.
Files containing the underlying binary data are placed into the `data_assets` directory.

## Reusability

This workflow is relatively tailored to the data which it was originally written
to analyse. Additional ensembles may be added to the analysis by adding relevant
files to the `raw_data` directory, and adding corresponding entries to the files
in the `metadata` directory. However, extending the analysis in this way has not
been as fully tested as the rest of the workflow, and is not guaranteed to be 
trivial for someone not already familiar with the code.

[datarelease]: https://doi.org/10.5281/zenodo.20111459
[miniforge]: https://github.com/conda-forge/miniforge
[paper]: https://doi.org/10.48550/arXiv.YYMM.XXXXX
[snakemake]: https://snakemake.github.io
[texlive]: https://tug.org/texlive/
