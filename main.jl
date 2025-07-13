using Pkg; Pkg.activate("src/src_jl"); Pkg.instantiate()

t0     = 0
deriv  = true
gevp   = false
use3x3 = true
plotting = true
symmetrise = false
average_equivalent_momenta = true

raw_path  = "./raw_data/"
plotpath  = "./assets/plots/"
tablepath = "./assets/tables/"
datapath  = "./data_assets/"

h5file_raw  = joinpath(datapath,"isospin1_sorted.hdf5")
h5file_com  = joinpath(datapath,"isospin1_merged.hdf5")
h5file_cor  = joinpath(datapath,"isospin1_corr.hdf5")
h5file_eig  = joinpath(datapath,"isospin1_eigenvalues.hdf5")
h5file_fit  = joinpath(datapath,"isospin1_fitresults.hdf5")
h5file_scat = joinpath(datapath,"isospin1_scattering.hdf5")
h5file_scat_fit = joinpath(datapath,"isospin1_fit_scatter.hdf5")

inputfiles = "metadata/input_files.csv"
infvolfile = "metadata/infinite_volume.csv"
fitparam   = "metadata/pipi_fitintervals.csv"
input_scatter = "metadata/scattering_input.csv"
input_scatter_fit = "metadata/fit_scatter_input.csv"
ensembles_list = "metadata/ensembles.csv"

overview_table    = joinpath(tablepath,"all_runs.csv")
analysed_table    = joinpath(tablepath,"analysed_runs.csv")

include("src/scripts_julia/plot_correlation_matrix.jl")
include("src/scripts_julia/plot_effective_masses.jl")
include("src/scripts_julia/plot_eigenvalues.jl")

run(`julia src/scripts_julia/parse_all_files.jl --path $(raw_path) --h5file $(h5file_raw) --inputfiles $(inputfiles)`)
run(`julia src/scripts_julia/write_tables.jl --h5file $(h5file_raw) --outfile $(overview_table)`)
run(`julia src/scripts_julia/write_tables.jl --h5file $(h5file_raw) --outfile $(analysed_table) --ensembles_list $(ensembles_list)`)
run(`julia src/scripts_julia/combine_runs.jl --h5file_in $(h5file_raw) --h5file_out $(h5file_com)`)
run(`julia src/scripts_julia/write_correlation_matrix.jl --h5file_in $(h5file_com) --h5file_out $(h5file_cor) --ensembles_list $(ensembles_list)`)
run(`julia src/scripts_julia/write_eigenvalues.jl --h5file_in $(h5file_cor) --h5file_out $(h5file_eig) --gevp $gevp --t0 $t0 --deriv $deriv --avg $average_equivalent_momenta --symmetrise $symmetrise`)
run(`python3 src/src_py/fitting.py $(h5file_eig) $(h5file_fit) $(fitparam)`)

plot_correlation_matrices(h5file_com,plotpath)
plot_eigenvalues(h5file_eig,plotpath)
plot_effective_masses(h5file_eig, h5file_fit, infvolfile, plotpath; use3x3)

ispath("tmp") || mkpath("tmp")
redirect_stdio(stdout="tmp/make.log",stderr="tmp/make.log") do 
    run(`bash src/zeta/compile.sh`)
end

cp(h5file_fit,h5file_scat,force=true)
run(`python3 src/src_py/scattering.py $(input_scatter) $(h5file_fit) $(h5file_scat)`)
cp(h5file_scat,h5file_scat_fit,force=true)
run(`python3 src/src_py/fit_scatter.py $(h5file_scat) $(h5file_scat_fit)`) 
run(`python3 src/src_py/plotting.py $(joinpath(plotpath,"scattering")) $(h5file_scat) $(h5file_scat_fit)`) 