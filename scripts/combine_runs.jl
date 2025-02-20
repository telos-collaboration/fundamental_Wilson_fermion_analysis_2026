using Pkg; Pkg.activate(".")
using HDF5
function check_matching_runs(file,name,key)
    fid  = h5open(file)[name]
    runs = keys(fid)
    vals = read.(Ref(fid),joinpath.(runs,key))
    @assert allequal(vals)
end
function check_lattice_params(file,name)
    for l in ["beta", "gauge group", "lattice", "quarkmasses", "Nconf", "configurations", "plaquette"]
        check_matching_runs(file,name,l)
    end
end
function copy_lattice_params(h5file_in,h5file_out,name)
    fid = h5open(h5file_in)[name]
    run_label = first(keys(fid))
    for l in ["beta", "gauge group", "lattice", "quarkmasses", "Nconf", "configurations", "plaquette"]
        h5write(h5file_out,joinpath(name,l),read(fid,joinpath(run_label,l)))
    end
end
function update_key(key)
    m = match(r"rho_g(?<g1>[1-3])(?<g2>[1-3])",key)
    return isnothing(m) ? key : "rho_g$(m.captures[1])_g$(m.captures[2])"
end
function old_key(key)
    m = match(r"rho_g(?<g1>[1-3])_g(?<g2>[1-3])",key)
    return isnothing(m) ? key : "rho_g$(m.captures[1])$(m.captures[2])"
end
function merge_runs(h5file,name)
    check_lattice_params(h5file,name)
end
function merge_runs(h5file_in, h5file_out, name )
    isfile(h5file_out) && rm(h5file_out)
    check_lattice_params(h5file_in,name)
    copy_lattice_params(h5file_in,h5file_out,name)

    fid  = h5open(h5file_in)[name]
    runs = keys(fid)

    p_ext = read.(Ref(fid),joinpath.(runs,"p_external"))
    p_ext_unique = unique(vcat(p_ext...))

    for p in p_ext_unique

        rids   = getindex.(Ref(fid),joinpath.(runs))
        labels = keys.(getindex.(rids,p))
        unique_labels = unique(vcat(labels...))

        for l in unique_labels

            # TODO: Check that the internal momenta always match
            p_internal = keys(first(rids)[joinpath(p,l)])

            for p_i in p_internal
                tmp_re = nothing 
                tmp_im = nothing 
                nsrc   = 0
                for r in rids
                    # The if-statement takes care of the labels for the off-diagonal
                    # vector meson correlators, since we have changed their label when
                    # we introduced the γiγ0 correlators 
                    if l ∈ keys(r[p]) 
                        re = read(r,joinpath(p,l,p_i,"C_re"))
                        im = read(r,joinpath(p,l,p_i,"C_im"))
                        tmp_re = isnothing(tmp_re) ? re : cat(tmp_re,re,dims=2)
                        tmp_im = isnothing(tmp_im) ? re : cat(tmp_im,im,dims=2)
                        nsrc += read(r,"Nsrc")
                    elseif old_key(l) in keys(r[p])
                        re = read(r,joinpath(p,old_key(l),p_i,"C_re"))
                        im = read(r,joinpath(p,old_key(l),p_i,"C_im"))  
                        tmp_re = isnothing(tmp_re) ? re : cat(tmp_re,re,dims=2)
                        tmp_im = isnothing(tmp_im) ? re : cat(tmp_im,im,dims=2)
                        nsrc += read(r,"Nsrc")
                    end
                end
                h5write(h5file_out,joinpath(name,p,l,p_i,"C_re"),tmp_re)
                h5write(h5file_out,joinpath(name,p,l,p_i,"C_im"),tmp_im)
                h5write(h5file_out,joinpath(name,p,l,p_i,"Nsrc"),nsrc)
            end
        end
    end
    h5write(h5file_out,joinpath(name,"p_external"),p_ext_unique)
end
h5file_in  = "data/isospin1_v2.hdf5"
h5file_out = "data/isospin1_merged.hdf5"
name   = "Lt24Ls14beta7.05m1-0.85m2-0.85"
merge_runs(h5file_in, h5file_out, name )
