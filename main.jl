using Pkg; Pkg.activate("."); Pkg.instantiate()
using ScatteringI1
using DelimitedFiles
using HDF5
using ProgressMeter
using LatticeUtils
using Plots
using LaTeXStrings
using PDFmerger
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

t0    = 8
deriv = true
path  = "/home/fabian/Dokumente/Physics/Data/DataVSC/measurements/"
plotpath = "./output/plots/"
datapath = "./output/data/"

h5file_raw = joinpath(datapath,"isospin1_sorted.hdf5")
h5file_com = joinpath(datapath,"isospin1_merged.hdf5")
h5file_cor = joinpath(datapath,"isospin1_corr.hdf5")
h5file_eig = joinpath(datapath,"isospin1_eigenvalues_t0_$(t0)_deriv_$deriv.hdf5")
h5file_fit = joinpath(datapath,"isospin1_fitresults_t0_$(t0)_deriv_$deriv.hdf5")

inputfiles = "input/input_files.csv"
infvolfile = "input/infinite_volume.csv"
fitparam   = "input/pipi_fitintervals.csv"

include("scripts/parse_all_files.jl")
include("scripts/combine_runs.jl")
include("scripts/write_correlation_matrix.jl")
include("scripts/write_eigenvalues.jl")
include("scripts/variational_analysis_meff.jl")

#parse_all_file(path,h5file_raw,inputfiles;single_file = true)
merge_all_runs(h5file_raw, h5file_com)
#write_correlation_matrix(h5file_com,h5file_cor)
#write_all_eigenvalues(h5file_cor,h5file_eig; t0, deriv, plotpath)
#run(`python3 scripts/fitting.py $(h5file_eig) $(h5file_fit) $(fitparam)`)
#plot_effective_masses(h5file_cor, h5file_fit, infvolfile, plotpath, fitparam; t0, deriv)