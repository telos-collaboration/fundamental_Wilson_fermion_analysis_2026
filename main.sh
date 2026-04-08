#!/bin/bash
set -e
export GKSwstype="100"

average_equivalent_momenta=true
num_resample_scattering=500

raw_path="./raw_data/"
plotpath="./assets/plots/"
plotpathscatter="./assets/plots/scattering/"
datapath="./data_assets/"
tablepath="./assets/tables/"

h5file_raw="data_assets/isospin1_sorted.hdf5"
h5file_com="data_assets/isospin1_merged.hdf5"
h5file_cor="data_assets/isospin1_corr.hdf5"
h5file_eig="data_assets/isospin1_eigenvalues.hdf5"
h5file_fit="data_assets/isospin1_fitresults.hdf5"
h5file_fit_full="data_assets/isospin1_fitresults_full.hdf5"
h5file_eig_gvp="data_assets/isospin1_eigenvalues_evp_gevp.hdf5"
h5file_fit_evp="data_assets/isospin1_fitresults_evp_gevp.hdf5"
# h5file_scat="data_assets/literature_scattering.hdf5"
h5file_scat="data_assets/isospin1_scattering.hdf5"
h5file_scat_fit="data_assets/isospin1_fit_scatter.hdf5"

inputfiles="metadata/input_files.csv"
infvolfile="metadata/infinite_volume.csv"
fitparam="metadata/pipi_fitintervals.csv"
fitparam_evp="metadata/pipi_fitintervals_evp_gevp.csv"
fitparam_meson="metadata/meson_fitintervals.csv"

input_scatter="metadata/scattering_input.csv"
fit_scatter_input="metadata/fit_scatter_input.csv"
ensembles_list="metadata/ensembles.csv"

# update all submodules in libs/
julia src/scripts_julia/instantiate.jl
julia src/scripts_julia/parse_all_files.jl --path $raw_path --h5file $h5file_raw --inputfiles $inputfiles
julia src/scripts_julia/write_tables.jl --h5file $h5file_raw --outfile "$tablepath/all_runs.csv"
julia src/scripts_julia/write_tables.jl --h5file $h5file_raw --outfile "$tablepath/analysed_runs.csv" --ensembles_list $ensembles_list
julia src/scripts_julia/combine_runs.jl --h5file_in $h5file_raw --h5file_out $h5file_com
julia src/scripts_julia/write_correlation_matrix.jl --h5file_in $h5file_com --h5file_out $h5file_cor --ensembles_list $ensembles_list
julia src/scripts_julia/write_eigenvalues.jl --h5file_in $h5file_cor --h5file_out $h5file_eig --metadata $fitparam --avg $average_equivalent_momenta
julia src/scripts_julia/write_eigenvalues.jl --h5file_in $h5file_cor --h5file_out $h5file_eig_gvp --metadata $fitparam_evp --avg $average_equivalent_momenta --swap_metadata "metadata/swaps.csv"
python3 src/src_py/fitting_eigenvalues.py $h5file_eig $h5file_fit $fitparam
cp $h5file_fit $h5file_fit_full
python3 src/src_py/fitting_mesons.py $h5file_eig $h5file_fit_full $fitparam_meson
python3 src/src_py/fitting_eigenvalues.py $h5file_eig_gvp $h5file_fit_evp $fitparam_evp
julia src/scripts_julia/write_table_fitresults.jl --h5file $h5file_fit_full --outfile "$tablepath/fit_results_3x3_tuned.csv"

julia src/scripts_julia/plot_diagrams.jl --h5file_in $h5file_com --plotpath $plotpath
julia src/scripts_julia/plot_eigenvalues.jl --h5file_in $h5file_eig --metadata $fitparam --plotname "$plotpath/eigenvalues.pdf"
julia src/scripts_julia/plot_eigenvalues.jl --h5file_in $h5file_eig_gvp --metadata $fitparam_evp --plotname "$plotpath/eigenvalues_evp_gevp.pdf"
julia src/scripts_julia/plot_eigenvalues_with_fits.jl --h5file_in $h5file_eig --plotpath $plotpath --metadata $fitparam --fitresults $h5file_fit_full
julia src/scripts_julia/plot_correlation_matrix_elements.jl --h5file_in $h5file_eig --plotpath $plotpath
julia src/scripts_julia/plot_meson_correlators.jl --h5file_in $h5file_eig --plotpath $plotpath --fitresults $h5file_fit_full
julia src/scripts_julia/plot_effective_masses.jl --h5file_eig $h5file_eig --h5file_fit $h5file_fit_full --plotpath $plotpath --infinite_volume $infvolfile --metadata $fitparam --plot_mesons true  --plotbasename "effective_masses"
julia src/scripts_julia/plot_effective_masses.jl --h5file_eig $h5file_eig_gvp --h5file_fit $h5file_fit_evp --plotpath $plotpath --infinite_volume $infvolfile --metadata $fitparam_evp --plot_mesons false --plotbasename "effective_masses_comparison"

mkdir -p tmp
bash libs/zeta/compile.sh  &> tmp/make.log

cp $h5file_fit_full $h5file_scat
python3 src/src_py/scattering.py $input_scatter $h5file_fit_full $h5file_scat $num_resample_scattering gauss
python3 src/src_py/plot_fit_scatter.py $plotpath/scattering $h5file_scat False

cp $h5file_scat $h5file_scat_fit
python3 src/src_py/fit_scatter.py $h5file_scat_fit $fit_scatter_input
python3 src/src_py/plot_fit_scatter.py $plotpath/scattering $h5file_scat_fit True

python3 src/src_py/E_L_plot.py $plotpath/scattering $h5file_scat_fit
python3 src/src_py/phase_shift_plots.py $plotpath/scattering $h5file_scat_fit
python3 src/src_py/cross_section_plot.py $plotpath/scattering $h5file_scat_fit

python3 src/src_py/result_tables.py $tablepath $h5file_scat_fit