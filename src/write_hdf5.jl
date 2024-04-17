function _write_lattice_setup(file,h5file;h5group="")
    h5write(h5file,joinpath(h5group,"plaquette"),plaquettes(file))
    h5write(h5file,joinpath(h5group,"configurations"),confignames(file))
    h5write(h5file,joinpath(h5group,"beta"),inverse_coupling(file))
    h5write(h5file,joinpath(h5group,"lattice"),latticesize(file))
    h5write(h5file,joinpath(h5group,"mass"),fermionmass(file))
end
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
function isospin1_to_hdf5(file,h5file;setup=true,ensemble="")
    setup &&  _write_lattice_setup(file,h5file;h5group=ensemble)
    Re, Im = parse_isospin_one(file)
    labels = label_list(file)
    
    Nlabels, Nconf, Nsrc, Np, T = size(Re)
    @assert length(labels) == Nlabels
    
    for i in 1:Nlabels
        channel, P_tot =_splitlabel(labels[i])
        pindex = channel=="d" ? (1:3) : (2:3)
        for px in pindex, py in pindex, pz in pindex
            p_diag = "p_diag($(px-2),$(py-2),$(pz-2))"
            h5label_re = joinpath(ensemble,P_tot,channel,p_diag,"C_re")
            h5label_im = joinpath(ensemble,P_tot,channel,p_diag,"C_im")
            
            h5write(h5file,h5label_re,Re[i,:,:,px,py,pz,:])
            h5write(h5file,h5label_im,Im[i,:,:,px,py,pz,:])
        end
    end
end