using DelimitedFiles: readdlm

function parse_all_files(path,h5file,inputfiles)
    println("Parse correlator data from raw log:")
    info  = readdlm(inputfiles,',',skipstart=1)
    for (name,file,run) in eachrow(info)

        file = joinpath(path,file)
        ens  = joinpath(name,run)

        dir = dirname(h5file)
        ispath(dir) || mkpath(dir) 
        isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true,sort=true,deduplicate=true)
    end
end