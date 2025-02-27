using Pkg; Pkg.activate(".",io=devnull)
using ScatteringI1
using DelimitedFiles
include("scripts/parse_all_files.jl")

h5file = "data/isospin1_sorted.hdf5"
path   = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"
path   = "/home/fabian/Dokumente/Physics/Data/DataVSC/measurements/"
info   = readdlm("input/input_files.csv",',',skipstart=1)
parse_all_file(path,h5file,info;single_file = true)

using Pkg; Pkg.activate(".",io=devnull)
using HDF5
include("scripts/combine_runs.jl")

h5file_in  = "data/isospin1_sorted.hdf5"
h5file_out = "data/isospin1_merged.hdf5"
ensembles  = keys(h5open(h5file_in))

isfile(h5file_out) && rm(h5file_out)
for ensemble in ensembles
    try 
        merge_runs(h5file_in, h5file_out, ensemble )
    catch
        @warn "Ensemble $ensemble cannot be merged"
        continue
    end
end

using Pkg; Pkg.activate(".",io=devnull)
using ScatteringI1
using ProgressMeter
using HDF5
include("scripts/write_correlation_matrix.jl")

file_in  = "data/isospin1_merged.hdf5"
file_out = "data/isospin1_corr.hdf5"
write_correlation_matrix(file_in,file_out)

#file_in2 = "data/isospin1_sorted.hdf5"
#file_out2= "data/isospin1_corr_allruns.hdf5"
#write_correlation_matrix(file_in2,file_out2;combined=false)

using Pkg; Pkg.activate(".",io=devnull)
using ScatteringI1
using HDF5
using LatticeUtils
using Plots
using LaTeXStrings
using ProgressMeter
using PDFmerger
include("scripts/write_eigenvalues.jl")
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

t0      = 8
deriv   = true
plotpath="./plots/"
infile  = "data/isospin1_corr.hdf5"
outfile = ifelse(deriv,"data/isospin1_eigenvalues_t0_$(t0)_deriv.hdf5","data/isospin1_eigenvalues_t0_$t0.hdf5")
write_all_eigenvalues(infile,outfile; t0, deriv)

run(`python3 scripts/fitting.py`)

using Pkg; Pkg.activate(".",io=devnull)
using HDF5
using Plots
using ScatteringI1
using LaTeXStrings
using LatticeUtils
using DelimitedFiles
using ProgressMeter
using PDFmerger
include("scripts/variational_analysis_meff.jl")
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

corr_file  = "data/isospin1_corr.hdf5"
fitresults = "data/isospin1_fitresults_t0_8_deriv.hdf5"
plotpath   = "plots/"
infvolfile = "input/infinite_volume.csv"
fitparam   = "input/pipi_fitintervals.csv"

deriv = true
t0    = 8

plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, fitparam; t0, deriv)