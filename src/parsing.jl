function _parse_data!(array,string;n=6) 
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    io = IOBuffer(string)
    for i in 1:n
        array[i] = Parsers.parse(Float64, io, opts)
    end
end
function _parse_channel_name(string)
    # split name at '_src_'. Everything before labels the measurement 
    # The following number labels the stochastic source
    pos = findfirst("_src_",string)
    f,l = first(pos), last(pos)
    n   = findnext('_',string,l+1)

    label = string[1:first(f)-1]
    src   = parse(Int,string[l+1:n-1])

    return label, src
end
function _count_labels(file)
    return length(label_list(file))
end
function label_list(file)
    cut  = length("[IO][0]")
    conf = 0
    labels = String[]
    for line in eachline(file)
        if startswith(line,"[IO][0]")
            startswith(line,"[IO][0]Configuration") && continue
            if isletter(line[cut+1])
                label, src = _parse_channel_name(line[cut+1:end])
                push!(labels,label) 
            end
        end
        # This overcounts by a factor of the number of sources
        if startswith(line,"[MAIN][0]Configuration from")
            conf += 1
            if conf == 2
                return unique(labels)
            end
        end
    end
    return unique(labels)
end
function _nconfs(file)
    nconf   = 0
    started = false 
    for line in eachline(file)
        if startswith(line,"[MAIN][0]Configuration from")
            started = true
        end
        if startswith(line,"[MAIN][0]Configuration : analysed")
            if started
                nconf += 1
                started = false
            else
                @error "Measurement not completed"
            end
        end
    end
    return nconf
end
function _sources(file)
    for line in eachline(file)
        if startswith(line,"[MAIN][0]num sources:")
            c = length("[MAIN][0]num sources:")
            n = findnext(',',line,c+1)
            Nsrc = parse(Int,line[c+1:n-1])
            return Nsrc 
        end
    end
end
function _find_pmax(file)
    p_max = typemin(Int)
    for line in eachline(file)
        if startswith(line,"[MAIN][0]The momenta are: ")
            for m in eachmatch(r"[0-9]+",line)
                p_max = max( parse(Int,m.match) , p_max)
            end
            @assert p_max > 0
            return p_max
        end
    end
end
function parse_isospin_one(file)
    T = first(latticesize(file))
    cut  = length("[IO][0]")
    # The measured  momenta are currently hard-coded in HiRep, i.e. the momenta are (0,0,0),(0,0,1),(0,1,1),(1,1,1) and permutations
    # plus negative momenta are allowed for the 4-point functions 
    pmax  = _find_pmax(file)
    Nmom  = 2pmax + 1 #(allowed values -pmax,...,0,...pmax) 
    Nconf = _nconfs(file)
    Nlab = _count_labels(file)
    Nsrc = _sources(file)
    tmp  = zeros(6) # temporary array for parsing result
    conf = 0 
    src  = 0
    li   = 0
    # fill arrays with NaNs. The idea is that not all momentum indices are used for all diagrams
    # All available entries will be replaced by finite Float64 numbers, the rest remains a NaN rather 
    # than a zero. 
    Re = fill(NaN,(Nlab,Nconf,Nsrc,Nmom,Nmom,Nmom,T))
    Im = fill(NaN,(Nlab,Nconf,Nsrc,Nmom,Nmom,Nmom,T))

    nlines = countlines(file)
    p = Progress(nlines)

    labels = label_list(file)
    for line in eachline(file)
        if startswith(line,"[IO][0]Configuration")
            if occursin("read",line)
                conf += 1
            end
            continue
        end
        if startswith(line,"[IO][0]")
            l = line[cut+1:end]
            # first line that starts here encodes the channel, source and configuration name
            if isletter(line[cut+1])
                label, src = _parse_channel_name(l)
                li = findfirst(isequal(label),labels)
            else
                _parse_data!(tmp, l)
                px, py, pz, t, re, im = tmp
                px, py, pz, t = Int(px), Int(py), Int(pz), Int(t)
                offset = pmax + 1
                # increase indices by pmax+1, to have one-based indexing
                # e.g. for pmax=1: index 1: p = -1
                #                  index 2: p =  0
                #                  index 3: p =  1
                Re[li,conf,src+1,px+offset,py+offset,pz+offset,t+1] = re
                Im[li,conf,src+1,px+offset,py+offset,pz+offset,t+1] = im
            end
        end
        next!(p)
    end
    finish!(p)
    return Re, Im
end