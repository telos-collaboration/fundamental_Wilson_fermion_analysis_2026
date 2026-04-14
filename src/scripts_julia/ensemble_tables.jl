using Pkg; Pkg.activate("src/src_jl")
using HDF5
using MadrasSokal
using Statistics

function ensemble_table(h5file)

    h5dset = h5open(h5file)
    ensembles = keys(h5dset)

    fields = 11
    table = Array{Any}(undef, (length(ensembles),fields))

    for (i,ens) in enumerate(ensembles)
        # get autocorrelation times
        P = read(h5dset[ens],"plaquette")
        Q = read(h5dset[ens],"Q")
        ind = read(h5dset[ens],"trajectory indices")
        Nskip = ind[2] - ind[1]
        τP = madras_sokal_time(P)
        τQ = madras_sokal_time(Q)

        # get lattice parameters
        Ncnf = length(read(h5dset[ens],"configurations"))
        beta = read(h5dset[ens],"beta")
        m0 = read(h5dset[ens],"quarkmass")[1]
        T,L = read(h5dset[ens],"lattice")[1:2]

        # ensemble type
        type = isapprox(m0,0.92) ? "heavy" : isapprox(m0,0.863) ? "medium" : "light" 
        table[i,1] = type
        table[i,2] = beta
        table[i,3] = m0
        table[i,4] = T
        table[i,5] = L
        table[i,6] = Ncnf
        table[i,7] = Nskip
        table[i,8] = MadrasSokal.errorstring(mean(P),std(P)/sqrt(Ncnf))
        table[i,9] = MadrasSokal.errorstring(mean(Q),std(Q)/sqrt(Ncnf))
        table[i,10] = MadrasSokal.errorstring(τP...)
        table[i,11] = MadrasSokal.errorstring(τQ...)
    end
    # sort table 
    table = sortslices(table,dims=1, by=x->x[[2,3,4]])
    # set up table footer and header
    header = raw"""
            \begin{tabular}{|l|c|c|c|c|c|c|c|c|c|c|}
                \hline \hline
                label & $\beta$ & $-am_0$ & $N_t$ & $N_s$ & $N_{\rm config.}$ & $N_{\rm skip}$ & $\langle P \rangle$ & $Q$ & $\tau^{\langle P \rangle}$ & $\tau^{Q}$ \\
            """
    footer = raw"""
                \hline \hline
            \end{tabular}
            """
    label = ""
    # actually write table to file            
    print(header)
    for r in eachrow(table)
        line = prod(string.(r) .* " & ")
        # remove last two characters to remove the trailing '&' and add a line brake
        line = line[1:end-2]*"\\\\"
        # if the label of the ensemble has changed insert extra hlines
        if label != r[1]
            println("\t\\hline\\hline")
        end
        label = r[1]
        println("\t"*line)
    end
    print(footer)
end

h5file = "data_assets/topology.hdf5"
ensemble_table(h5file)