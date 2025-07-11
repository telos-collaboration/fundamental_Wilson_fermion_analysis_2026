function parse_all_file(path,h5file,inputfiles;single_file = true)
    println("Parse correlator data from raw log:")
    info  = readdlm(inputfiles,',',skipstart=1)
    for (name,file,run) in eachrow(info)

        file = joinpath(path,file)
        ens  = joinpath(name,run)

        if single_file
            dir = dirname(h5file)
            ispath(dir) || mkpath(dir) 
            isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true,sort=true,deduplicate=true)
        else
            dir = joinpath(dirname(h5file),"ensembles")
            ispath(dir) || mkpath(dir) 
            f = joinpath(dir,"$(name)_$run.hdf5")
            isospin1_to_hdf5(file,f;ensemble="",setup=true,sort=true,deduplicate=true)
        end
    end
end