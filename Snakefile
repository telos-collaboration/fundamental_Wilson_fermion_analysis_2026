t0=0
deriv="true"
gevp="false"
use3x3="true"
plotting="true"
symmetrise="true"
average_equivalent_momenta="true"

plotpath="assets/plots"
datapath="data_assets"
tablepath="assets/tables"

h5file_fit="data_assets/isospin1_fitresults.hdf5"
h5file_scat="data_assets/isospin1_scattering.hdf5"
h5file_scat_fit="data_assets/isospin1_fit_scatter.hdf5"

infvolfile="metadata/infinite_volume.csv"
input_scatter="metadata/scattering_input.csv"
input_scatter_fit="metadata/fit_scatter_input.csv"

rule all:
    input: 
        h5file="data_assets/isospin1_fitresults.hdf5",
        fittable=f"{tablepath}/fit_results_3x3.csv",
        plots=[
            "assets/plots/diagrams.pdf",
            "assets/plots/diagrams_3x3.pdf",
            "assets/plots/eigenvalues.pdf",
            "assets/plots/effective_masses_(g)evp.pdf",
        ],
        plots_sigma=expand("{p}/sigma1_{ens}_fit_True.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"]),
        plots_cot=expand("{p}/p3cotPS_{ens}_fit_True.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"]),
        plots_EL=expand("{p}/E_L_{ens}_levels_{b}.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"], b=[True, False]),
        plots_ECML=expand("{p}/E_CM_L_{ens}_levels_{b}.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"], b=[True, False]),
        binary="libs/zeta/out/get_wlm.so",
        scatterh5="data_assets/isospin1_scattering.hdf5",
        scatterfith5="data_assets/isospin1_fit_scatter.hdf5",

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
        table="{tablepath}/fit_results_3x3.csv"
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --outfile {output.table}'

rule plot_correlation_matrices:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/plot_correlation_matrix.jl",
        h5file="data_assets/isospin1_merged.hdf5",
    output: 
        plots=[
            "assets/plots/diagrams.pdf",
            "assets/plots/diagrams_3x3.pdf",
            ]
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --plotpath {plotpath}'

rule plot_eigenvalues:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/plot_eigenvalues.jl",
        h5file="data_assets/isospin1_eigenvalues.hdf5",
    output: 
        plots=[
            "assets/plots/eigenvalues.pdf",
            ]
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --plotpath {plotpath} --three_by_three {use3x3}'

rule plot_correlation_matrix_elements:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/plot_correlation_matrix_elements.jl",
        h5file="data_assets/isospin1_eigenvalues.hdf5",
    output: 
        plots=[
            "assets/plots/correlation_matrix_elements.pdf",
            ]
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --plotpath {plotpath}'

rule plot_effective_masses:
    input:
        julia_instantiated="tmp/julia_ready",
        script="src/scripts_julia/plot_effective_masses.jl",
        h5file_eig="data_assets/isospin1_eigenvalues.hdf5",
        h5file_fit="data_assets/isospin1_fitresults.hdf5",
        infinite_volume="metadata/infinite_volume.csv",
    output:
        plots=[
            "assets/plots/effective_masses_(g)evp.pdf",
            ]
    conda:
        "environment.yml"
    shell:
        'julia {input.script} --h5file_eig {input.h5file_eig} --h5file_fit {input.h5file_fit} --infinite_volume {input.infinite_volume} --plotpath {plotpath} --three_by_three {use3x3}'


rule compile_zeta:
    input:
        code="libs/zeta/"
    output:
        binary="libs/zeta/out/get_wlm.so"
    conda:
        "environment.yml"
    shell:
        'bash libs/zeta/compile.sh  &> tmp/make.log'

rule scattering_analysis:
    input:
        h5file="data_assets/isospin1_fitresults.hdf5",
        input_scatter="metadata/scattering_input.csv",
        script="src/src_py/scattering.py",
    output:
        h5file="data_assets/isospin1_scattering.hdf5",
    conda:
        "environment.yml"
    shell: """ 
        cp {input.h5file} {output.h5file}
        python3 {input.script} {input.input_scatter} {input.h5file} {output.h5file}
        """ 

rule scattering_fits:
    input:
        h5file="data_assets/isospin1_scattering.hdf5",
        script="src/src_py/fit_scatter.py",
    output:
        h5file="data_assets/isospin1_fit_scatter.hdf5",
    conda:
        "environment.yml"
    shell: """ 
        cp {input.h5file} {output.h5file}
        python3 {input.script} {input.h5file} {output.h5file}
        """ 

rule scattering_fit_plots:
    input:
        h5file_scatter="data_assets/isospin1_scattering.hdf5",
        h5file_scatter_fit="data_assets/isospin1_fit_scatter.hdf5",
        script="src/src_py/plotting.py",
    output:
        plots_sigma=expand("{p}/sigma1_{ens}_fit_True.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"]),
        plots_cot=expand("{p}/p3cotPS_{ens}_fit_True.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"]),
        plots_EL=expand("{p}/E_L_{ens}_levels_{b}.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"], b=[True, False]),
        plots_ECML=expand("{p}/E_CM_L_{ens}_levels_{b}.pdf", p=f"{plotpath}/scattering", ens=["res","close_res","non_res"], b=[True, False]),
    conda:
        "environment.yml"
    shell: """ 
        python3 {input.script} {plotpath}/scattering {input.h5file_scatter} {input.h5file_scatter_fit}
        """ 
