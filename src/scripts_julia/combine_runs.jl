using HDF5: h5open, h5write

function check_matching_runs(file,ensemble,key;verbose=true)
    fid  = h5open(file)[ensemble]
    runs = keys(fid)
    vals = read.(Ref(fid),joinpath.(runs,key))
    if !allequal(vals) && verbose
        @warn "mismatch in $ensemble for $key: $vals"
    end
    @assert allequal(vals) "mismatch in $ensemble for $key"
end
function check_lattice_params(file,ensemble)
    for l in ["beta", "gauge group", "lattice", "quarkmasses", "Nconf", "configurations", "plaquette"]
        check_matching_runs(file,ensemble,l)
    end
end
function copy_lattice_params(h5file_in,h5file_out,ensemble)
    fid = h5open(h5file_in)[ensemble]
    run_label = first(keys(fid))
    for l in ["beta", "gauge group", "lattice", "quarkmasses", "Nconf", "configurations", "plaquette"]
        h5write(h5file_out,joinpath(ensemble,l),read(fid,joinpath(run_label,l)))
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
# This obatins any internal momenta, external momenta or diagram labels present in the hdf5 file
function unique_labels_momenta(file_id)
    runs = keys(file_id)
    # Arrays that store all labels & momenta encountered
    all_p_internal = AbstractString[]
    all_labels = AbstractString[] 
    # First get all external momenta
    p_external = read.(Ref(file_id),joinpath.(runs,"p_external"))
    p_external_unique = unique(vcat(p_external...))
    # Then find all momenta & labels 
    for r in runs
        for p_e in p_external_unique
            l1 = joinpath(r,p_e)
            if haskey(file_id,l1)
                labels = keys(file_id[l1])
                append!(all_labels,labels)
                for lab in labels
                    l2 = joinpath(r,p_e,lab)
                    p_internal = keys(file_id[l2])
                    append!(all_p_internal,p_internal)
                end
            end
        end
    end
    return p_external_unique, unique(update_key.(all_labels)), unique(all_p_internal)
end
function merge_runs(h5file_in, h5file_out, ensemble )
    check_lattice_params(h5file_in,ensemble)
    copy_lattice_params(h5file_in,h5file_out,ensemble)

    fid  = h5open(h5file_in)[ensemble]
    runs = keys(fid)
    rids = getindex.(Ref(fid),joinpath.(runs))
    p_ext_unique, unique_labels, p_internal_unique = unique_labels_momenta(fid)

    for p in p_ext_unique
        for l in unique_labels
            for p_i in p_internal_unique
                tmp_re = nothing 
                tmp_im = nothing 
                nsrc   = 0
                for r in rids
                    # The if-statement takes care of the labels for the off-diagonal
                    # vector meson correlators, since we have changed their label when
                    # we introduced the γiγ0 correlators 
                    if haskey(r,p) && haskey(r[p],l)
                        if  p_i in keys(r[joinpath(p,l)])
                            re = read(r,joinpath(p,l,p_i,"C_re"))
                            im = read(r,joinpath(p,l,p_i,"C_im"))
                            tmp_re = isnothing(tmp_re) ? re : cat(tmp_re,re,dims=2)
                            tmp_im = isnothing(tmp_im) ? im : cat(tmp_im,im,dims=2)
                            nsrc += read(r,"Nsrc")
                        end
                    elseif haskey(r,p) && haskey(r[p],old_key(l)) 
                        if  p_i in keys(r[joinpath(p,old_key(l))])
                            re = read(r,joinpath(p,old_key(l),p_i,"C_re"))
                            im = read(r,joinpath(p,old_key(l),p_i,"C_im"))  
                            tmp_re = isnothing(tmp_re) ? re : cat(tmp_re,re,dims=2)
                            tmp_im = isnothing(tmp_im) ? im : cat(tmp_im,im,dims=2)
                            nsrc += read(r,"Nsrc")
                        end
                    end
                end
                if !isnothing(tmp_re) && !isnothing(tmp_im) && nsrc > 0
                    h5write(h5file_out,joinpath(ensemble,p,l,p_i,"C_re"),tmp_re)
                    h5write(h5file_out,joinpath(ensemble,p,l,p_i,"C_im"),tmp_im)
                    h5write(h5file_out,joinpath(ensemble,p,l,p_i,"Nsrc"),nsrc)
                end
            end
        end
    end
    h5write(h5file_out,joinpath(ensemble,"p_external"),p_ext_unique)
end
function merge_all_runs(h5file_in, h5file_out)
    isfile(h5file_out) && rm(h5file_out)
    ensembles  = keys(h5open(h5file_in))
    for ensemble in ensembles
        try 
            merge_runs(h5file_in, h5file_out, ensemble )
        catch
            @warn "Ensemble $ensemble cannot be merged"
            continue
        end
    end
end
