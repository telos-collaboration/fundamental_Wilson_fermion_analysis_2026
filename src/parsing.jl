function latticesize(file)
    for line in eachline(file)
        if occursin("Global size is",line)
            pos  = last(findfirst("Global size is",line))+1
            sizestring  = lstrip(line[pos:end])
            latticesize = parse.(Int,split(sizestring,"x"))
            return latticesize
        end
    end
end
function _parse_data!(array,string) 
    substrings = eachsplit(string," ",keepempty=false)
    for (i,s) in enumerate(substrings)
        array[i] = parse(Float64,s)
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
    return length(_label_list(file))
end
function _label_list(file)
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
function _parse_isospin_one(file)
    T = first(latticesize(file))
    cut  = length("[IO][0]")
    # The measured  momenta are currently hard-coded in HiRep, i.e. the momenta are (0,0,0),(0,0,1),(0,1,1),(1,1,1) and permutations
    # plus negative momenta are allowed for the 4-point functions 
    Nmom  = 3 #(allowed values -1,0,1) 
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
    Re = zeros(Nlab,Nconf,Nsrc,Nmom,Nmom,Nmom,T) .* NaN
    Im = zeros(Nlab,Nconf,Nsrc,Nmom,Nmom,Nmom,T) .* NaN

    labels = _label_list(file)
    for line in eachline(file)
        if startswith(line,"[IO][0]")
            if startswith(line,"[IO][0]Configuration")
                continue
            end
            l = line[cut+1:end]
            # first line that starts here encodes the channel, source and configuration name
            if isletter(line[cut+1])
                label, src = _parse_channel_name(l)
                li = findfirst(isequal(label),labels)
            else
                _parse_data!(tmp, l)
                px, py, pz, t, re, im = tmp
                px, py, pz, t = Int(px), Int(py), Int(pz), Int(t)
                # increase indices by one, to have one-based indexing
                # for momenta: index 1: p = -1
                #              index 2: p =  0
                #              index 3: p =  1
                Re[li,conf,src+1,px+2,py+2,pz+2,t+1] = re
                Im[li,conf,src+1,px+2,py+2,pz+2,t+1] = im
            end
        end
        if startswith(line,"[MAIN][0]Configuration from")
            conf += 1
        end
    end
    return Re, Im
end
