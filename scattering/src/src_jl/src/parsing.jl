function _parse_reim(io,opts) 
    re = Parsers.parse(Float64,io, opts)
    im = Parsers.parse(Float64,io, opts)
    return re, im
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
    io = HiRepParsing.makestream(file)
    for line in eachline(io)
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
                close(io)
                return unique(labels)
            end
        end
    end
    close(io)
    return unique(labels)
end
function _nconfs(file)
    nconf   = 0
    started = false
    io = HiRepParsing.makestream(file)
    for line in eachline(io)
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
    close(io)
    return nconf
end
function _sources(file)
    io = HiRepParsing.makestream(file)
    for line in eachline(io)
        if startswith(line,"[MAIN][0]num sources:")
            c = length("[MAIN][0]num sources:")
            n = findnext(',',line,c+1)
            Nsrc = parse(Int,line[c+1:n-1])
            close(io)
            return Nsrc 
        end
    end
    close(io)
end
function _mom_from_label(label)
    label == "pi" && (return [0,0,0])
    if contains(label,"p0") 
        return [0,0,0]
    else 
        return _parse_momentum(label)
    end
end
parse_isospin_one(file;kws...) = parse_isospin_one(file,_nconfs(file);desc="Progress:")
function parse_isospin_one(file,Nconf;desc="Progress:")
    T = first(latticesize(file))
    cut  = length("[IO][0]")
    # The measured  momenta are currently hard-coded in HiRep, i.e. the momenta are (0,0,0),(0,0,1),(0,1,1),(1,1,1) and permutations
    # plus negative momenta are allowed for the 4-point functions 
    Nlab = _count_labels(file)
    Nsrc = _sources(file)
    conf = 0 
    src  = 0

    # fill arrays with NaNs. The idea is that not all momentum indices are used for all diagrams
    # All available entries will be replaced by finite Float64 numbers, the rest remains a NaN rather 
    # than a zero. 
    Re = fill(NaN,(Nlab,Nconf,Nsrc,2,T))
    Im = fill(NaN,(Nlab,Nconf,Nsrc,2,T))
    
    # store the current label, its index among all labels and the associated external momentum
    li    = 0
    label = ""
    p_ext = [0,0,0]

    # set up options for Parsers
    opts = Parsers.Options(delim=' ', ignorerepeated=true)
    p = Progress(Nconf; desc)

    labels = label_list(file)
    io = HiRepParsing.makestream(file)
    for line in eachline(io)
        if startswith(line,"[IO][0]Configuration")
            if occursin("read",line)
                conf += 1
                next!(p)
            end
            continue
        end
        if startswith(line,"[IO][0]")
            l = SubString(line,cut+1)
            # first line that starts here encodes the channel, source and configuration name
            if isletter(line[cut+1])
                label, src = _parse_channel_name(l)
                li = findfirst(isequal(label),labels)
                p_ext = _mom_from_label(label)
            else
                # create IO buffer from current line
                io = IOBuffer(l)
                # Parsing the following integers is the most expensive part of the entire parsing step 
                # Here, I check after every momentum component if we can continue before parsing additional momentum components
                px = Parsers.parse(Int64, io, opts)
                px != p_ext[1] && px != 0 && continue
                py = Parsers.parse(Int64, io, opts)
                py != p_ext[2] && py != 0 && continue
                pz = Parsers.parse(Int64, io, opts)
                pz != p_ext[3] && pz != 0 && continue
                # (px,py,pz) == (0,0,0) save in index (1)
                # (px,py,pz) == p_ext   save in index (2)
                if px==0 && py==0 && pz==0
                    t = Parsers.parse(Int64,io,opts)
                    re, im = _parse_reim(io,opts) 
                    Re[li,conf,src+1,1,t+1] = re
                    Im[li,conf,src+1,1,t+1] = im
                elseif px == p_ext[1] && py == p_ext[2] && pz == p_ext[3] 
                    t = Parsers.parse(Int64,io,opts)
                    re, im = _parse_reim(io,opts) 
                    Re[li,conf,src+1,2,t+1] = re
                    Im[li,conf,src+1,2,t+1] = im
                end
            end
        end
    end
    finish!(p)
    close(io)
    return Re, Im
end