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
function fermionmasses(file;pattern="[MAIN][0]Mass[0]")
    masses = Float64[]
    for line in eachline(file)
        if occursin(pattern,line)
            s = split(line,(","))
            for i in eachindex(s)
                m = parse(Float64,split(s[i],"=")[2])
                append!(masses,m)
            end
            return masses
        end
    end
end
function plaquettes(file)
    plaquettes = Float64[]
    for line in eachline(file)
        if occursin("Plaquette",line)
            line = replace(line,"="=>" ")
            line = replace(line,":"=>" ")
            p = parse(Float64,split(line)[end])
            append!(plaquettes,p)
        end
    end
    return plaquettes
end
function _match_config_name(filename)
    regex = r".*/(?<run>[^/]*)_(?<T>[0-9]+)x(?<L>[0-9]+)x[0-9]+x[0-9]+nc[0-9]+(?:r[A-Z]+)?(?:nf[0-9]+)?b(?<beta>[0-9]+\.[0-9]+)?(?:m-?[0-9]+\.[0-9]+)?n(?<conf>[0-9]+)"
    return match(regex,filename)
end
function inverse_coupling(file)
    try
        l = split(file,"beta")[end]
        β = parse(Float64,split(l,"m")[1])
        return β
    catch
        for line in eachline(file)
            if occursin("Configuration from",line)
                match = _match_config_name(line)
                β = parse(Float64,match[:beta])
                return β
            end
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
# increase indices by one, to have one-based indexing
# encode the momenta in a linear index. This is essentially base 3
#      p_index = (px+1)(py+1)(pz+1)_base3 
_momentum_to_index(px,py,pz) = 9*(px+1) + 3*(py+1) + (pz+1) + 1
_index_to_momentum(index) = reverse(Tuple(digits(index-1,base=3,pad=3).-1))
_index_to_momentum_label(index) = string(_index_to_momentum(index))
function parse_isospin_one(file)
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
    #
    # Furthermore, I encode the momena in base three, since onlc three possible values are hardcoded
    Re = zeros(Nlab,Nconf,Nsrc,Nmom^3,T) .* NaN
    Im = zeros(Nlab,Nconf,Nsrc,Nmom^3,T) .* NaN

    labels = label_list(file)
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
                p_index = _momentum_to_index(px,py,pz)
                Re[li,conf,src+1,p_index,t+1] = re
                Im[li,conf,src+1,p_index,t+1] = im
            end
        end
        if startswith(line,"[MAIN][0]Configuration from")
            conf += 1
        end
    end
    return Re, Im
end
function confignames(file)
    fns = AbstractString[]
    for line in eachline(file)
        if occursin("read",line)
            if occursin("Configuration",line)
                pos1 = findlast('/',line)
                pos2 = findnext(']',line,pos1)
                push!(fns,line[pos1+1:pos2-1])
            end
        end
    end
    return fns
end