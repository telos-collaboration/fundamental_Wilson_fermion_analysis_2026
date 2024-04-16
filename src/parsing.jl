_parse_channel_name(x) = nothing

function _parse_data!(array,string) 
    substrings = eachsplit(string," ",keepempty=false)
    for (i,s) in enumerate(substrings)
        array[i] = parse(Float64,s)
    end
end

function parse_isospin_one(file)
    cut  = length("[IO][0]")
    tmp  = zeros(6) # temporary array for parsing result
    conf = 0 
    
    for line in eachline(file)
        if startswith(line,"[IO][0]")
            if startswith(line,"[IO][0]Configuration")
                continue
            end
            l = line[cut+1:end]
            # first line that starts here encode the channel, source and configuration name
            if isletter(line[cut+1])
                @show l
                _parse_channel_name(l)
                break
            else
                _parse_data!(tmp, l)
                break
            end
        end
        if startswith(line,"[MAIN][0]Configuration from")
            conf += 1
            @show conf
        end
    end

end

