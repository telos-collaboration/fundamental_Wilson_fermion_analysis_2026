function _splitlabel(label)
    label == "pi" && (return label, "p(0,0,0)")
    if contains(label,"p0") 
        l = replace(label,"_p0"=>"")
        return l, "p(0,0,0)"
    elseif contains(label,"p(")
        f = first(findfirst("_p(",label))
        n = findlast(')',label)
        l = label[1:f-1]*label[n+1:end]
        p = label[f+1:n]
        return l, p
    end
end
unique_indices(v) = unique(i -> v[i], eachindex(v))
function isospin1_to_hdf5(file,h5file;setup=true,ensemble="",sort=true,deduplicate=true)
    setup &&  HiRepParsing._write_lattice_setup(file,h5file;h5group=ensemble,smearing=false,sort,deduplicate)
    perm   = h5read(h5file,joinpath(ensemble,"sort_permutation"))
    inds   = h5read(h5file,joinpath(ensemble,"deduplicated_indices"))
    Nconf0 = length(perm)
    Nconf  = length(inds)

    Re, Im = parse_isospin_one(file,Nconf0;desc=ensemble)
    Nlabels, Nconf0, Nsrc, Nmom, T = size(Re)
    labels = label_list(file)
    @assert Nmom == 2 
    @assert length(labels) == Nlabels
    
    p_external = unique(last.(_splitlabel.(labels)))
    h5write(h5file,joinpath(ensemble,"Nsrc"),Nsrc)
    h5write(h5file,joinpath(ensemble,"Nconf"),Nconf)
    h5write(h5file,joinpath(ensemble,"p_external"),p_external)

    tmpRe = zeros(Nconf,Nsrc,T)
    tmpIm = zeros(Nconf,Nsrc,T)
    for i in 1:Nlabels

        channel, P_tot =_splitlabel(labels[i])
        # (px,py,pz) == (0,0,0) have been saved in index (1)
        # (px,py,pz) == p_ext   have been saved in index (2)
        h5label_re = joinpath(ensemble,P_tot,channel,"p_diag(0,0,0)","C_re")
        h5label_im = joinpath(ensemble,P_tot,channel,"p_diag(0,0,0)","C_im")            
        @. tmpRe = Re[i,perm[inds],:,1,:]
        @. tmpIm = Im[i,perm[inds],:,1,:]
        h5write(h5file,h5label_re,tmpRe)
        h5write(h5file,h5label_im,tmpIm)

        if P_tot != "p(0,0,0)"
            h5label_re = joinpath(ensemble,P_tot,channel,"p_diag$(P_tot[2:end])","C_re")
            h5label_im = joinpath(ensemble,P_tot,channel,"p_diag$(P_tot[2:end])","C_im")            
            @. tmpRe = Re[i,perm[inds],:,2,:]
            @. tmpIm = Im[i,perm[inds],:,2,:]
            h5write(h5file,h5label_re,tmpRe)
            h5write(h5file,h5label_im,tmpIm)
        end

    end
end