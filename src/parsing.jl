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

    # increase number of sources by one, to have one-based indexing
    return label, src+1
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
function parse_isospin_one(file)
    cut  = length("[IO][0]")
    Nmom = 2 # This is currently hard-coded in HiRep, i.e. the momenta are (0,0,0),(0,0,1),(0,1,1),(1,1,1) and permutations
    Nlab = _count_labels(file)
    Nsrc = _sources(file)
    tmp  = zeros(6) # temporary array for parsing result
    tmpRe = zeros(T,Nmom,Nmom,Nmom,Nsrc,Nlab)
    tmpIm = zeros(T,Nmom,Nmom,Nmom,Nsrc,Nlab)
    conf = 0 
    for line in eachline(file)
        if startswith(line,"[IO][0]")
            if startswith(line,"[IO][0]Configuration")
                continue
            end
            l = line[cut+1:end]
            # first line that starts here encodes the channel, source and configuration name
            if isletter(line[cut+1])
                label, src = _parse_channel_name(l)
                @show label, src 
            else
                _parse_data!(tmp, l)
                px, py, pz, t, re, im = tmp
            end
        end
        if startswith(line,"[MAIN][0]Configuration from")
            conf += 1
            @show conf
        end
    end
end
