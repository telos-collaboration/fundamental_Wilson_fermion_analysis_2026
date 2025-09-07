#!/bin/bash
set -e

average_equivalent_momenta=true
num_resample_scattering=200

version="yannick"

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
h5file_fit_evp="data_assets/isospin1_fitresults_evp.hdf5"
# h5file_scat="data_assets/literature_scattering.hdf5"
h5file_scat="data_assets/isospin1_scattering.hdf5"
h5file_scat_fit="data_assets/isospin1_fit_scatter.hdf5"

inputfiles="metadata/input_files.csv"
infvolfile="metadata/infinite_volume.csv"
fitparam="metadata/pipi_fitintervals.csv"
fitparam_evp="metadata/pipi_fitintervals_gevp_test.csv"
fitparam_meson="metadata/meson_fitintervals.csv"

input_scatter="metadata/scattering_input.csv"
fit_scatter_input="metadata/fit_scatter_input.csv"
ensembles_list="metadata/ensembles.csv"

# update all submodules in libs/
git submodule update --init --recursive
julia src/scripts_julia/instantiate.jl
julia src/scripts_julia/parse_all_files.jl --path $raw_path --h5file $h5file_raw --inputfiles $inputfiles
julia src/scripts_julia/write_tables.jl --h5file $h5file_raw --outfile "$tablepath/all_runs.csv"
julia src/scripts_julia/write_tables.jl --h5file $h5file_raw --outfile "$tablepath/analysed_runs.csv" --ensembles_list $ensembles_list
julia src/scripts_julia/combine_runs.jl --h5file_in $h5file_raw --h5file_out $h5file_com
julia src/scripts_julia/write_correlation_matrix.jl --h5file_in $h5file_com --h5file_out $h5file_cor --ensembles_list $ensembles_list
julia src/scripts_julia/write_eigenvalues.jl --h5file_in $h5file_cor --h5file_out $h5file_eig --metadata $fitparam --avg $average_equivalent_momenta
python3 src/src_py/fitting_eigenvalues.py $h5file_eig $h5file_fit $fitparam
python3 src/src_py/fitting_mesons.py $h5file_eig $h5file_fit $fitparam_meson
python3 src/src_py/fitting_eigenvalues.py $h5file_eig $h5file_fit_evp $fitparam_evp
julia src/scripts_julia/write_table_fitresults.jl --h5file $h5file_fit --outfile "$tablepath/fit_results_3x3_tuned.csv"

julia src/scripts_julia/plot_diagrams.jl --h5file_in $h5file_com --plotpath $plotpath
julia src/scripts_julia/plot_eigenvalues.jl --h5file_in $h5file_eig --plotpath $plotpath --metadata $fitparam
julia src/scripts_julia/plot_eigenvalues_with_fits.jl --h5file_in $h5file_eig --plotpath $plotpath --metadata $fitparam --fitresults $h5file_fit
julia src/scripts_julia/plot_correlation_matrix_elements.jl --h5file_in $h5file_eig --plotpath $plotpath
julia src/scripts_julia/plot_meson_correlators.jl --h5file_in $h5file_eig --plotpath $plotpath --fitresults $h5file_fit
julia src/scripts_julia/plot_effective_masses.jl --h5file_eig $h5file_eig --h5file_fit $h5file_fit --plotpath $plotpath --infinite_volume $infvolfile --metadata $fitparam

mkdir -p tmp
bash libs/zeta/compile.sh  &> tmp/make.log

cp $h5file_fit $h5file_scat
python3 src/src_py/scattering.py $input_scatter $h5file_fit $h5file_scat $num_resample_scattering gauss
python3 src/src_py/plot_scatter.py $plotpathscatter $h5file_scat $version
cp $h5file_scat $h5file_scat_fit
#python3 src/src_py/fit_scatter.py $h5file_scat $h5file_scat_fit $fit_scatter_input
#python3 src/src_py/plotting.py $plotpath/scattering $h5file_scat_fit