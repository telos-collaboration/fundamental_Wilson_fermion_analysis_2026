function _parse_data!(array::Array{T},io;n) where T 
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    for i in 1:n
        array[i] = Parsers.parse(T, io, opts)
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
function _mom_from_label(label)
    label == "pi" && (return [0,0,0])
    if contains(label,"p0") 
        return [0,0,0]
    else 
        return _parse_momentum(label)
    end
end
function parse_isospin_one(file;desc="Progress:")
    T = first(latticesize(file))
    cut  = length("[IO][0]")
    # The measured  momenta are currently hard-coded in HiRep, i.e. the momenta are (0,0,0),(0,0,1),(0,1,1),(1,1,1) and permutations
    # plus negative momenta are allowed for the 4-point functions 
    Nconf = _nconfs(file)
    Nlab = _count_labels(file)
    Nsrc = _sources(file)
    conf = 0 
    src  = 0
    tmpInt = zeros(Int64,4) # temporary array for parsing result
    tmpFlt = zeros(2) # temporary array for parsing result

    # fill arrays with NaNs. The idea is that not all momentum indices are used for all diagrams
    # All available entries will be replaced by finite Float64 numbers, the rest remains a NaN rather 
    # than a zero. 
    Re = fill(NaN,(Nlab,Nconf,Nsrc,2,T))
    Im = fill(NaN,(Nlab,Nconf,Nsrc,2,T))
    
    # store the current label, its index among all labels and the associated external momentum
    li    = 0
    label = ""
    p_ext = [0,0,0]

    nlines = countlines(file)
    p = Progress(nlines; desc)

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
                p_ext = _mom_from_label(label)
            else
                io = IOBuffer(l) # create IO buffer from current line
                _parse_data!(tmpInt,io;n=4) 
                px, py, pz, t = tmpInt
                if any(isnothing,tmpInt) || isnothing(src) || isnothing(li)
                    @error "line could not be parsed correctly" line label li src conf
                end  
                # (px,py,pz) == (0,0,0) save in index (1)
                # (px,py,pz) == p_ext   save in index (2)
                if (px,py,pz) == (0,0,0)
                    _parse_data!(tmpFlt,io;n=2)
                    re, im = tmpFlt 
                    Re[li,conf,src+1,1,t+1] = re
                    Im[li,conf,src+1,1,t+1] = im
                elseif [px,py,pz] == p_ext
                    _parse_data!(tmpFlt,io;n=2)
                    re, im = tmpFlt 
                    Re[li,conf,src+1,2,t+1] = re
                    Im[li,conf,src+1,2,t+1] = im
                end
            end
        end
        next!(p)
    end
    finish!(p)
    return Re, Im
end