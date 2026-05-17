import os
os.environ["GKSwstype"] = "100"


rule julia_instantiate:
    input:
        script="scattering/src/scripts_julia/instantiate.jl",
    output:
        julia_instantiated="intermediary_data/julia_ready",
    conda: 
        "../envs/scattering.yml"
    shell:
        "julia {input.script} && touch {output.julia_instantiated}"


rule parse_logs:
    input:
        script="scattering/src/scripts_julia/parse_all_files.jl",
        julia_instantiated="intermediary_data/julia_ready",
        metadata="metadata/scattering/input_files.csv",
    output:
        h5file="intermediary_data/scattering/isospin1_sorted.hdf5",
    conda:
        "../envs/scattering.yml"
    # Start parsing early,
    # since it is time consuming and many other processes depend on it
    priority: 10
    shell:
        'julia {input.script} --path raw_data/scattering/ --h5file {output.h5file} --inputfiles {input.metadata}'


rule ensembles_tables:
    input:
        script="scattering/src/scripts_julia/write_tables.jl",
        julia_instantiated="intermediary_data/julia_ready",
        h5file="intermediary_data/scattering/isospin1_sorted.hdf5",
        metadata="metadata/scattering/ensembles.csv"
    output:
        table="assets/scattering/tables/analysed_runs.csv"
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --ensembles_list {input.metadata} --outfile {output.table}'


rule combine_runs_hdf5:
    input:
        script="scattering/src/scripts_julia/combine_runs.jl",
        julia_instantiated="intermediary_data/julia_ready",
        h5file_in="intermediary_data/scattering/isospin1_sorted.hdf5",
    output:
        h5file_out="intermediary_data/scattering/isospin1_merged.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --h5file_out {output.h5file_out}'


rule correlation_matrix:
    input:
        script="scattering/src/scripts_julia/write_correlation_matrix.jl",
        julia_instantiated="intermediary_data/julia_ready",
        h5file_in="intermediary_data/scattering/isospin1_merged.hdf5",
        metadata="metadata/scattering/ensembles.csv"
    output:
        h5file_out="data_assets/scattering/isospin1_corr.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --h5file_out {output.h5file_out} --ensembles_list {input.metadata}'


rule write_eigenvalues:
    input:
        script="scattering/src/scripts_julia/write_eigenvalues.jl",
        julia_instantiated="intermediary_data/julia_ready",
        h5file_in="data_assets/scattering/isospin1_corr.hdf5",
        metadata="metadata/scattering/pipi_fitintervals.csv",
    output:
        h5file_out="data_assets/scattering/isospin1_eigenvalues.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --h5file_out {output.h5file_out}  --metadata {input.metadata} --avg true'


rule write_eigenvalues_evp_vs_gevp:
    input:
        script="scattering/src/scripts_julia/write_eigenvalues.jl",
        julia_instantiated="intermediary_data/julia_ready",
        h5file_in="data_assets/scattering/isospin1_corr.hdf5",
        metadata="metadata/scattering/pipi_fitintervals_evp_gevp.csv",
        swaps="metadata/scattering/swaps.csv",
    output:
        h5file_out="data_assets/scattering/isospin1_eigenvalues_evp_gevp.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --h5file_out {output.h5file_out}  --metadata {input.metadata}  --swap_metadata {input.swaps} --avg true'


rule fit_eigenvalues:
    input:
        script="scattering/src/src_py/fitting_eigenvalues.py",
        metadata="metadata/scattering/pipi_fitintervals.csv",
        h5file_in="data_assets/scattering/isospin1_eigenvalues.hdf5",
    output:
        h5file_out="intermediary_data/scattering/isospin1_fitresults.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'python3 {input.script} {input.h5file_in} {output.h5file_out} {input.metadata}'


rule fit_eigenvalues_comparison:
    input:
        script="scattering/src/src_py/fitting_eigenvalues.py",
        metadata="metadata/scattering/pipi_fitintervals_evp_gevp.csv",
        h5file_in="data_assets/scattering/isospin1_eigenvalues_evp_gevp.hdf5",
    output:
        h5file_out="data_assets/scattering/isospin1_fitresults_evp_gevp.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        'python3 {input.script} {input.h5file_in} {output.h5file_out} {input.metadata}'


rule fit_mesons:
    input:
        script="scattering/src/src_py/fitting_mesons.py",
        metadata="metadata/scattering/meson_fitintervals.csv",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues.hdf5",
        h5file_fit_pipi="intermediary_data/scattering/isospin1_fitresults.hdf5",
    output:
        h5file_out="intermediary_data/scattering/isospin1_full_fitresults.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        cp {input.h5file_fit_pipi} {output.h5file_out}
        python3 {input.script} {input.h5file_eig} {output.h5file_out} {input.metadata}
        '''


rule fit_result_table:
    input:
        script="scattering/src/scripts_julia/write_table_fitresults.jl",
        h5file="intermediary_data/scattering/isospin1_full_fitresults.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        table="assets/scattering/tables/fit_results_3x3_tuned.csv"
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file {input.h5file} --outfile {output.table}'


rule plot_individual_diagrams:
    input:
        script="scattering/src/scripts_julia/plot_diagrams.jl",
        h5file="intermediary_data/scattering/isospin1_merged.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plots2x2="assets/scattering/plots/diagrams.pdf",
        plots3x3="assets/scattering/plots/diagrams_3x3.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --plotpath "assets/scattering/plots/"'


rule plot_eigenvalues:
    input:
        script="scattering/src/scripts_julia/plot_eigenvalues.jl",
        h5file="data_assets/scattering/isospin1_eigenvalues.hdf5",
        metadata="metadata/scattering/pipi_fitintervals.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot="assets/scattering/plots/eigenvalues.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --metadata {input.metadata} --plotname {output.plot}'


rule plot_eigenvalues_comparison:
    input:
        script="scattering/src/scripts_julia/plot_eigenvalues.jl",
        h5file="data_assets/scattering/isospin1_eigenvalues_evp_gevp.hdf5",
        metadata="metadata/scattering/pipi_fitintervals_evp_gevp.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot="assets/scattering/plots/eigenvalues_evp_gevp.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file} --metadata {input.metadata} --plotname {output.plot}'

rule plot_eigenvalues_with_fits:
    input:
        script="scattering/src/scripts_julia/plot_eigenvalues_with_fits.jl",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues.hdf5",
        h5file_fit="intermediary_data/scattering/isospin1_fitresults.hdf5",
        metadata="metadata/scattering/pipi_fitintervals.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot="assets/scattering/plots/eigenvalues_with_fit.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_eig} --plotpath "assets/scattering/plots/" --metadata {input.metadata} --fitresults {input.h5file_fit}'


rule plot_correlation_matrix_elements:
    input:
        script="scattering/src/scripts_julia/plot_correlation_matrix_elements.jl",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot="assets/scattering/plots/correlation_matrix_elements.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_eig} --plotpath "assets/scattering/plots/"'


rule plot_meson_correlators:
    input:
        script="scattering/src/scripts_julia/plot_meson_correlators.jl",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues.hdf5",
        h5file_fit="intermediary_data/scattering/isospin1_fitresults.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot="assets/scattering/plots/meson_correlators.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_eig} --plotpath "assets/scattering/plots/" --fitresults {input.h5file_fit}'


rule plot_effective_masses:
    input:
        script="scattering/src/scripts_julia/plot_effective_masses.jl",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues.hdf5",
        h5file_fit="intermediary_data/scattering/isospin1_full_fitresults.hdf5",
        fitparam="metadata/scattering/pipi_fitintervals.csv",
        infinite_volume="metadata/scattering/infinite_volume.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot_gevp="assets/scattering/plots/effective_masses_(g)evp.pdf",
        plot_mesons="assets/scattering/plots/effective_masses_mesons.pdf",
        plot_mesons_p0="assets/scattering/plots/effective_masses_mesons_p0.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_eig {input.h5file_eig} --h5file_fit {input.h5file_fit} --infinite_volume {input.infinite_volume} --metadata {input.fitparam} --plotpath "assets/scattering/plots" --plot_mesons true --plotbasename "effective_masses"'


rule plot_effective_masses_comparison:
    input:
        script="scattering/src/scripts_julia/plot_effective_masses.jl",
        h5file_eig="data_assets/scattering/isospin1_eigenvalues_evp_gevp.hdf5",
        h5file_fit="data_assets/scattering/isospin1_fitresults_evp_gevp.hdf5",
        fitparam="metadata/scattering/pipi_fitintervals_evp_gevp.csv",
        infinite_volume="metadata/scattering/infinite_volume.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        plot_gevp="assets/scattering/plots/effective_masses_comparison_(g)evp.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_eig {input.h5file_eig} --h5file_fit {input.h5file_fit} --infinite_volume {input.infinite_volume} --metadata {input.fitparam} --plotpath "assets/scattering/plots" --plot_mesons false --plotbasename "effective_masses_comparison"'


rule compile_zeta:
    input:
        script="libs/zeta/compile.sh",
    output:
        binary="libs/zeta/out/get_wlm.so",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        mkdir -p intermediary_data/
        bash {input.script}  &> intermediary_data/make.log
        '''


rule luescher_analysis:
    input:
        script="scattering/src/src_py/scattering.py",
        h5file_fit="intermediary_data/scattering/isospin1_full_fitresults.hdf5",
        metadata="metadata/scattering/scattering_input.csv",
        binary="libs/zeta/out/get_wlm.so",
    output:
        h5file_out="intermediary_data/scattering/isospin1_scattering.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        cp {input.h5file_fit} {output.h5file_out}
        python3 {input.script} {input.metadata} {input.h5file_fit} {output.h5file_out} 500 gauss
        '''


rule plot_luescher_results:
    input:
        script="scattering/src/src_py/plot_fit_scatter.py",
        h5file="intermediary_data/scattering/isospin1_scattering.hdf5",
    output:
        plot1 = "assets/scattering/plots/scattering/p3cotPS_vs_p2star_heavy.pdf",
        plot2 = "assets/scattering/plots/scattering/p3cotPS_vs_p2star_medium.pdf",
        plot3 = "assets/scattering/plots/scattering/p3cotPS_Ecm_vs_s_medium.pdf",
        plot4 = "assets/scattering/plots/scattering/p3cotPS_Ecm_vs_s_light.pdf",
        plot5 = "assets/scattering/plots/scattering/PS_heavy.pdf",
        plot6 = "assets/scattering/plots/scattering/PS_medium.pdf",
        plot7 = "assets/scattering/plots/scattering/PS_light.pdf",
        plot8 = "assets/scattering/plots/scattering/sigma_heavy.pdf",
        plot9 = "assets/scattering/plots/scattering/sigma_medium.pdf",
        plot10 = "assets/scattering/plots/scattering/sigma_light.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/plots/scattering {input.h5file} False
        '''


rule plot_luescher_fits:
    input:
        script="scattering/src/src_py/plot_fit_scatter.py",
        h5file="data_assets/scattering/isospin1_fit_scatter.hdf5",
    output:
        plot1 = "assets/scattering/plots/scattering/p3cotPS_vs_p2star_heavy_ERE0.pdf",
        plot2 = "assets/scattering/plots/scattering/p3cotPS_vs_p2star_heavy_ERE1.pdf",
        plot3 = "assets/scattering/plots/scattering/p3cotPS_Ecm_vs_s_light_BWI.pdf",
        plot4 = "assets/scattering/plots/scattering/p3cotPS_Ecm_vs_s_light_BWII.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/plots/scattering {input.h5file} True
        '''


rule fit_phase_shifts:
    input:
        script="scattering/src/src_py/fit_scatter.py",
        h5file_in="intermediary_data/scattering/isospin1_scattering.hdf5",
        metadata="metadata/scattering/fit_scatter_input.csv",
    output:
        h5file_out="data_assets/scattering/isospin1_fit_scatter.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        cp {input.h5file_in} {output.h5file_out}
        python3 {input.script} {output.h5file_out} {input.metadata}
        '''


rule volume_dependence_plots:
    input:
        script="scattering/src/src_py/E_L_plot.py",
        h5file="intermediary_data/scattering/isospin1_scattering.hdf5",
    output:
        plot_heavy = "assets/scattering/plots/scattering/E_vs_L_heavy.pdf",
        plot_lighter = "assets/scattering/plots/scattering/E_vs_L_medium_light.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/plots/scattering {input.h5file} 
        '''


rule phase_shift_plots:
    input:
        script="scattering/src/src_py/phase_shift_plots.py",
        h5file="data_assets/scattering/isospin1_fit_scatter.hdf5",
    output:
        plot_light = "assets/scattering/plots/scattering/phase_shift_plot_light.pdf",
        plot_medium = "assets/scattering/plots/scattering/phase_shift_plot_medium.pdf",
        plot_heavy = "assets/scattering/plots/scattering/phase_shift_plot_heavy.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/plots/scattering {input.h5file} 
        '''


rule crossection_plots:
    input:
        script="scattering/src/src_py/cross_section_plot.py",
        h5file="data_assets/scattering/isospin1_fit_scatter.hdf5",
        I2results="raw_data/14_dim/scattering_Fig5.3_b6.900_m-0.920.hdf5",
    output:
        cross_section_plot2 = "assets/scattering/plots/scattering/sigma_comb_units.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/plots/scattering {input.h5file} 
        '''


rule luescher_tables:
    input:
        script="scattering/src/src_py/result_tables.py",
        h5file="data_assets/scattering/isospin1_fit_scatter.hdf5",
    output:
        table_heavy = "assets/scattering/tables/latex_scattering_b6.90_m-0.920_table.txt", 
        table_medium = "assets/scattering/tables/latex_scattering_b7.05_m-0.863_table.txt", 
        table_light = "assets/scattering/tables/latex_scattering_b7.05_m-0.867_table.txt",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        python3 {input.script} assets/scattering/tables/ {input.h5file} 
        '''


rule parse_flows:
    input:
        script="scattering/src/src_py/package_flows.py",
    output:
        h5file="data_assets/scattering/topology.hdf5",
    conda:
        "../envs/scattering.yml"
    shell:
        '''
        for filename in raw_data/gradient_flow/*/topology/out/out_flow; do 
            name=$(echo $filename | grep -o "Lt.*FUN")
            echo $name 
            python3 {input.script} --h5_filename {output.h5file} --ensemble $name $filename
        done
        '''


rule plot_meson_finite_volume:
    input:
        script="scattering/src/scripts_julia/finite_volume.jl",
        h5file_in="intermediary_data/scattering/isospin1_full_fitresults.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        table = "intermediary_data/finite_volume.csv",
        fv_light = "assets/scattering/plots/finite_volume/fv_b7.05_m0.867.pdf",
        fv_medium = "assets/scattering/plots/finite_volume/fv_b7.05_m0.863.pdf",
        fv_heavy = "assets/scattering/plots/finite_volume/fv_b6.9_m0.92.pdf",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --table_out {output.table} --plotpath "assets/scattering/plots/finite_volume/" '


rule ensemble_table:
    input:
        script="scattering/src/scripts_julia/ensemble_tables.jl",
        h5file_in="data_assets/scattering/topology.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        table = "assets/scattering/tables/table_II.tex",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --table_out {output.table} '


rule phase_shift_csv:
    input:
        script="scattering/src/scripts_julia/phase_shift_table.jl",
        h5file_in="data_assets/scattering/isospin1_fit_scatter.hdf5",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        table = "intermediary_data/phase_shift_{pattern}.csv",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --h5file_in {input.h5file_in} --outfile {output.table} --pattern {wildcards.pattern} '


rule phase_shift_tex:
    input:
        script="scattering/src/scripts_julia/phase_shift_tex.jl",
        table = "intermediary_data/phase_shift_{pattern}.csv",
        metadata = "metadata/scattering/fit_scatter_input.csv",
        julia_instantiated="intermediary_data/julia_ready",
    output:
        table = "assets/scattering/tables/{pattern}.tex",
    conda:
        "../envs/scattering.yml"
    shell:
        'julia {input.script} --csv_in {input.table} --outfile {output.table} --metadata {input.metadata}'

