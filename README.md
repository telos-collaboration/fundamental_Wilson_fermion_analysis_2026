# Vector meson scattering in Sp(4) gauge theory—Analysis workflow

The workflow in this repository performs the analysis presented in the paper

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- LaTeX, for example, from [TeX Live][texlive]

## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/fzierler/ScatteringI1
   cd ScatteringI1
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

[miniforge]: https://github.com/conda-forge/miniforge
[texlive]: https://tug.org/texlive/
