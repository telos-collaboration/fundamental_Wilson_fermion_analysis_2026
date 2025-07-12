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

overview_table    = joinpath(tablepath,"all_runs.csv")
analysed_table    = joinpath(tablepath,"analysed_runs.csv")

include("src/scripts_julia/parse_all_files.jl")
include("src/scripts_julia/combine_runs.jl")
include("src/scripts_julia/write_correlation_matrix.jl")
include("src/scripts_julia/write_eigenvalues.jl")
include("src/scripts_julia/plot_correlation_matrix.jl")
include("src/scripts_julia/plot_effective_masses.jl")
include("src/scripts_julia/plot_eigenvalues.jl")
include("src/scripts_julia/write_tables.jl")

only_ens = nothing
only_ens = [
            "Lt24Ls14beta6.9m-0.92",
            "Lt32Ls16beta6.9m-0.92",
            "Lt32Ls24beta6.9m-0.92",
            "Lt36Ls16beta7.05m-0.863",
            "Lt36Ls24beta7.05m-0.863",
            "Lt36Ls36beta7.05m-0.863",
            "Lt36Ls16beta7.05m-0.867",
            "Lt36Ls24beta7.05m-0.867",
            "Lt36Ls36beta7.05m-0.867",
        ]

parse_all_files(raw_path,h5file_raw,inputfiles)
all_runs_table(h5file_raw,overview_table;)
all_runs_table(h5file_raw,analysed_table;only_ens)
merge_all_runs(h5file_raw, h5file_com)

write_correlation_matrix(h5file_com,h5file_cor;only_ens)
plot_correlation_matrices(h5file_com,plotpath;only_ens)
write_all_eigenvalues(h5file_cor,h5file_eig; t0, deriv, gevp, average_equivalent_momenta,symmetrise)
plot_eigenvalues(h5file_eig,plotpath)
run(`python3 src/src_py/fitting.py $(h5file_eig) $(h5file_fit) $(fitparam)`)
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