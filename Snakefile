t0=0
deriv="true"
gevp="false"
use3x3="true"
plotting="true"
symmetrise="true"
average_equivalent_momenta="true"

plotpath="./assets/plots/"
datapath="./data_assets/"

tablepath="assets/tables/"

h5file_fit="data_assets/isospin1_fitresults.hdf5"
h5file_scat="data_assets/isospin1_scattering.hdf5"
h5file_scat_fit="data_assets/isospin1_fit_scatter.hdf5"

infvolfile="metadata/infinite_volume.csv"
input_scatter="metadata/scattering_input.csv"
input_scatter_fit="metadata/fit_scatter_input.csv"

rule all:
    input: 
        h5file="data_assets/isospin1_fitresults.hdf5",
        fittable=f"{tablepath}/fit_results_3x3_tuned.csv",

rule julia_instantiate:
    input:
        script="src/scripts_julia/instantiate.jl",
    output:
        julia_instantiated="tmp/julia_ready",
    conda:
        "environment.yml"
    shell:
        """
        git submodule update --init --recursive
        julia {input.script}
        touch {output.julia_instantiated}
        """

rule parse_hdf5:
    input:
        julia_instantiated="tmp/julia_ready",
        metadata="metadata/input_files.csv",
        script="src/scripts_julia/parse_all_files.jl",
        path="raw_data/",
    output:
        h5file="data_assets/isospin1_sorted.hdf5",
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --path {input.path} --h5file {output.h5file} --inputfiles {input.metadata}'

rule ensemble_table:
    input:
        julia_instantiated="tmp/julia_ready",
        h5file="data_assets/isospin1_sorted.hdf5",
        script="src/scripts_julia/write_tables.jl",
    output:
        table="{tablepath}/all_runs.csv",
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {output.h5file} --outfile {output.table}'

rule restricted_ensemble_table:
    input:
        julia_instantiated="tmp/julia_ready",
        h5file="data_assets/isospin1_sorted.hdf5",
        script="src/scripts_julia/write_tables.jl",
        ensembles_list="metadata/ensembles.csv",
    output:
        table="{tablepath}/analysed_runs.csv",
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {output.h5file} --outfile {output.table} --ensemble_list {input.ensembles_list}'

rule combine_runs:
    input:
        julia_instantiated="tmp/julia_ready",
        h5file="data_assets/isospin1_sorted.hdf5",
        script="src/scripts_julia/combine_runs.jl",
    output:
        h5file="data_assets/isospin1_merged.hdf5"
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --h5file_out {output.h5file}'

rule write_correlation_matrices:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/write_correlation_matrix.jl",
        h5file="data_assets/isospin1_merged.hdf5",
        ensembles_list="metadata/ensembles.csv",
    output:
        h5file="data_assets/isospin1_corr.hdf5",
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --h5file_out {output.h5file} --ensembles_list {input.ensembles_list}'

rule write_eigenvalues:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/write_eigenvalues.jl",
        h5file="data_assets/isospin1_corr.hdf5",
        ensembles_list="metadata/ensembles.csv",
    output:
        h5file="data_assets/isospin1_eigenvalues.hdf5",
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --h5file_out {output.h5file} --gevp {gevp} --t0 {t0} --deriv {deriv} --avg {average_equivalent_momenta} --symmetrise {symmetrise}'

rule fit_eigenvalues:
    input:
        h5file="data_assets/isospin1_eigenvalues.hdf5",
        script="src/src_py/fitting.py",
        metadata="metadata/pipi_fitintervals_3x3.csv",
    output:
        h5file="data_assets/isospin1_fitresults.hdf5",
    conda:
        "environment.yml"
    shell:
        'python3 {input.script} {input.h5file} {output.h5file} {input.metadata}'

rule fitresult_table:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/write_table_fitresults.jl",
        h5file="data_assets/isospin1_fitresults.hdf5",
    output: 
        table="{tablepath}/fit_results_3x3_tuned.csv"
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --outfile {output.table}'

#julia src/scripts_julia/plot_correlation_matrix.jl --h5file_in $h5file_com --plotpath $plotpath
#julia src/scripts_julia/plot_eigenvalues.jl --h5file_in $h5file_eig --plotpath $plotpath --three_by_three $use3x3
#julia src/scripts_julia/plot_effective_masses.jl --h5file_eig $h5file_eig --h5file_fit $h5file_fit --plotpath $plotpath --infinite_volume $infvolfile --three_by_three $use3x3
#
#mkdir -p tmp
#bash libs/zeta/compile.sh  &> tmp/make.log
#
#cp $h5file_fit $h5file_scat
#python3 src/src_py/scattering.py $input_scatter $h5file_fit $h5file_scat
#cp $h5file_scat $h5file_scat_fit
#python3 src/src_py/fit_scatter.py $h5file_scat $h5file_scat_fit
#python3 src/src_py/plotting.py $plotpath/scattering $h5file_scat $h5file_scat_fit
