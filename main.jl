using Pkg; Pkg.activate("."); Pkg.instantiate()
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

deriv  = true
gevp   = false
t0     = 0
use3x3 = true
average_equivalent_momenta = false

path  = "/home/fabian/Documents/Physics/Data/"
path  = "/home/fabian/Dokumente/Physics/Data/"
plotpath  = "./output/plots/"
datapath  = "./output/data/"
scatpath  = "./output/scattering/"
tablepath = "./output/tables/"

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

inputfiles = "input/input_files.csv"
infvolfile = "input/infinite_volume.csv"
fitparam   = "input/pipi_fitintervals.csv"

overview_table    = joinpath(tablepath,"all_runs.csv")
analysed_table    = joinpath(tablepath,"analysed_runs.csv")
yannick_fmt_table = joinpath(tablepath,"yannick_format_t0_$(t0)_deriv_$deriv.dat")

include("scripts/parse_all_files.jl")
include("scripts/combine_runs.jl")
include("scripts/write_correlation_matrix.jl")
include("scripts/write_eigenvalues.jl")
include("scripts/variational_analysis_meff.jl")
include("scripts/write_tables.jl")

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
half_sources=false

parse_all_file(path,h5file_raw,inputfiles;single_file = true)
all_runs_table(h5file_raw,overview_table;)
all_runs_table(h5file_raw,analysed_table;only_ens)
merge_all_runs(h5file_raw, h5file_com)

write_correlation_matrix(h5file_com,h5file_cor;plotpath,plotting,only_ens)
write_all_eigenvalues(h5file_cor,h5file_eig; t0, deriv, plotpath, plotting, use3x3, gevp, average_equivalent_momenta)
run(`python3 scripts/fitting.py $(h5file_eig) $(h5file_fit) $(fitparam)`)
plot_effective_masses(h5file_cor, h5file_fit, infvolfile, plotpath, fitparam; t0, deriv, gevp, use3x3, half_sources, average_equivalent_momenta)

redirect_stdio(stdout="make.log",stderr="make.log") do 
    run(`bash rho_pipi_scattering_analysis/zeta/compile.sh`)
end

cd("rho_pipi_scattering_analysis")
cp("../$(h5file_fit)","../$(h5file_scat)",force=true)
run(`python3 src/scattering.py test`)
cp("../$(h5file_scat)","../$(h5file_scat_fit)",force=true)
run(`python3 src/fit_scatter.py`) 
run(`python3 src/plotting.py`) 