using Pkg; Pkg.activate("scattering/src/src_jl")
using HDF5: h5open, h5write
using ScatteringI1
using LatticeUtils
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",frame=:box,legend=:left,markeralpha=0.9,size=(400,300),titlefontsize=11)

function format_fit_names(fits)
    # get dimension of correlation matrix and t0 
    rx = r"(?<N>[0-9])x([0-9])|t0_(?<t0>[0-9]+)"
    fit_fmt = copy(fits)

    # collect N and t0 for later sorting
    param = zeros(Int,length(fits)) 

    for (i,f) in enumerate(fits)
        m = match(rx,f)
        N = m["N"]
        t0 = m["t0"]
        if isnothing(t0)
            fit_fmt[i] = L"evp $(%$N\times %$N)$"
            param[i] = typemax(Int)
        else
            fit_fmt[i] = L"$t_0=%$t0$"
            param[i] = parse(Int,t0)
        end
    end
    return fit_fmt, param
end

h5file = "data_assets/scattering/isospin1_fitresults_evp_gevp.hdf5"
h5id = h5open(h5file)
ensembles = keys(h5id)

for ens in ensembles

    # only look at datasets corresponding to fits
    fits = filter(contains("evp"),keys(h5id[ens]))
    moms = filter(!isequal("p(0,0,0)"),keys(h5id[ens][fits[1]]))
    fits_fmt, param = format_fit_names(fits)
    perm = sortperm(param)

    for p in moms
        E  = hcat([read(h5id, joinpath(ens,f,p,"A1/E")) for f in fits ]...)
        ΔE = hcat([read(h5id, joinpath(ens,f,p,"A1/Delta_E")) for f in fits ]...)

        # create a title
        rx = r"Lt(?<Lt>[0-9]+)Ls(?<Ls>[0-9]+)beta(?<beta>[0-9]+.[0-9]+)m(?<m>-[0-9]+.[0-9]+)"
        m = match(rx,ens)
        L, T, β, m0 = m["Lt"], m["Ls"], m["beta"], m["m"]
        p_fmt = replace(p,"p"=>"")
        title = L"{%$T}\times{%$L}^3, \beta=%$β, am_0={%$m0}, \mathbf{d}=%$(p_fmt)"

        # obtain tick label    

        plt = plot(;title)
        scatter!(plt,E[1,perm],yerr=ΔE[1,perm], label = "state #0")
        scatter!(plt,E[2,perm],yerr=ΔE[2,perm], label = "state #1")
        plot!(plt,xticks=(1:length(fits),fits_fmt[perm]),xrotation = 45)
        display(plt)
    end

end