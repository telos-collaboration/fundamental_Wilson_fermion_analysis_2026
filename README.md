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

3. For now, I host the input on both VSC4 and DiaL3. They are rsyncable zstd
archives. Download and decompress the zstd archive from either of the following two locations.

    ```shellsession
    /scratch/dp208/shared/raw_data.tar.zst
    ``` 
    ```shellsession
    /gpfs/data/fs71564/zierler/raw_data.tar.zst
    ``` 

4. Create and activate the conda environment.

    ```shellsession
    conda env create -f environment.yml
    conda activate rho-pi-pi
    ``` 
5. Run the workflow

    ```shellsession
    bash main.sh
    ```

On an AMD 5950X CPU the current workflow takes roughly 25 minutes to complete.

## Output

Output plots, tables, and definitions are placed into the `data_assets` directory.

## Reusability

This workflow is relatively tailored to the data which it was originally written
to analyse. Additional ensembles may be added to the analysis by adding relevant
files to the `raw_data` directory, and adding corresponding entries to the files
in the `metadata` directory. However, extending the analysis in this way has not
been as fully tested as the rest of the workflow, and is not guaranteed to be 
trivial for someone not already familiar with the code.

[miniforge]: https://github.com/conda-forge/miniforge
[texlive]: https://tug.org/texlive/
