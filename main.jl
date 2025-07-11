using Pkg; Pkg.activate("src/src_jl"); Pkg.instantiate()
using ScatteringI1
using DelimitedFiles
using HDF5
using ProgressMeter
using LatticeUtils
using Plots
using LaTeXStrings
using PDFmerger
using Statistics
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

t0     = 0
deriv  = true
gevp   = false
use3x3 = true
symmetrise = false
average_equivalent_momenta = true

path  = "./raw_data/"
plotpath  = "./data_assets/plots/"
datapath  = "./data_assets/data/"
scatpath  = "./data_assets/scattering/"
tablepath = "./data_assets/tables/"

h5file_raw = joinpath(datapath,"isospin1_sorted.hdf5")
h5file_com = joinpath(datapath,"isospin1_merged.hdf5")
h5file_cor = joinpath(datapath,"isospin1_corr.hdf5")
if gevp
    h5file_eig      = joinpath(datapath,"isospin1_eigenvalues_gevp_t0_$(t0)_deriv_$deriv.hdf5")
    h5file_fit      = joinpath(datapath,"isospin1_fitresults_gevp_t0_$(t0)_deriv_$deriv.hdf5")
    h5file_scat     = joinpath(scatpath,"isospin1_scattering_gevp_t0_$(t0)_deriv_$deriv.hdf5")
    h5file_scat_fit = joinpath(scatpath,"isospin1_fit_scatter_gevp_t0_$(t0)_deriv_$deriv.hdf5")
else
    h5file_eig      = joinpath(datapath,"isospin1_eigenvalues_evp_deriv_$deriv.hdf5")
    h5file_fit      = joinpath(datapath,"isospin1_fitresults_evp_deriv_$deriv.hdf5")
    h5file_scat     = joinpath(scatpath,"isospin1_scattering_evp_deriv_$deriv.hdf5")
    h5file_scat_fit = joinpath(scatpath,"isospin1_fit_scatter_evp_deriv_$deriv.hdf5")
end

inputfiles = "metadata/input_files.csv"
infvolfile = "metadata/infinite_volume.csv"
fitparam   = "metadata/pipi_fitintervals.csv"

overview_table    = joinpath(tablepath,"all_runs.csv")
analysed_table    = joinpath(tablepath,"analysed_runs.csv")
yannick_fmt_table = joinpath(tablepath,"yannick_format_t0_$(t0)_deriv_$deriv.dat")

include("src/scripts_julia/parse_all_files.jl")
include("src/scripts_julia/combine_runs.jl")
include("src/scripts_julia/write_correlation_matrix.jl")
include("src/scripts_julia/write_eigenvalues.jl")
include("src/scripts_julia/variational_analysis_meff.jl")
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

plotting = true

#parse_all_file(path,h5file_raw,inputfiles;single_file = true)
all_runs_table(h5file_raw,overview_table;)
all_runs_table(h5file_raw,analysed_table;only_ens)
merge_all_runs(h5file_raw, h5file_com)

write_correlation_matrix(h5file_com,h5file_cor;only_ens)
plot_correlation_matrices(h5file_com,plotpath;only_ens)
write_all_eigenvalues(h5file_cor,h5file_eig; t0, deriv, gevp, average_equivalent_momenta,symmetrise)
plot_eigenvalues(h5file_eig,plotpath)
run(`python3 src/src_py/fitting.py $(h5file_eig) $(h5file_fit) $(fitparam)`)
plot_effective_masses(h5file_cor, h5file_fit, infvolfile, plotpath, fitparam; t0, deriv, gevp, use3x3, average_equivalent_momenta, symmetrise)

ispath("tmp") || mkpath("tmp")
redirect_stdio(stdout="tmp/make.log",stderr="tmp/make.log") do 
    run(`bash src/zeta/compile.sh`)
end

ispath(scatpath) || mkpath(scatpath)
cd("src")
cp("../$(h5file_fit)","../$(h5file_scat)",force=true)
run(`python3 src_py/scattering.py test`)
cp("../$(h5file_scat)","../$(h5file_scat_fit)",force=true)
run(`python3 src_py/fit_scatter.py`) 
run(`python3 src_py/plotting.py`) 