using DelimitedFiles
using LaTeXStrings
include("utils.jl")

function row_colouring(x)
    x == "Graz" && (return "gray!25")
    x == "2019Singlets" && (return "blue!25")
    x == "2019NoSinglets" && (return "green!25")
    x == "Hasenbusch" && (return "yellow!25")
end
function write_overview(infile="input/overview.csv",outfile="output/tex_table.tex")
    data, header = readdlm(infile,';',header=true)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    io = open(outfile,"w")

    # write header
    table_layout  = repeat("|c",10)*"|"
    write(io,"\\begin{tabular}{$table_layout}\n\t\\hline\n")
    write(io,"\t",L"$\beta$ & $m^0_f$ & $L$ & $T$ & $N_{\rm conf}$ & $m_{PS}$ & $m_{V}$ & $m_{S}$ & $m_{AV}$  & $m_{\eta'}$ \\ \hline","\n")

    nrow, ncol = size(data)

    for (j,row) in enumerate(eachrow(data))
        # determine colour of the row based on origin of ensemble
        write(io,"\t\t\\rowcolor{$(row_colouring(last(row)))}\t")
        # the first five entries are exact parameters of the ensemble
        for i in 1:5
            write(io,"$(row[i]) & ")
        end
        # the remaining ones are masses, that come with uncertainties
        for i in 6:2:14
            if !isempty(row[i])
                write(io,"$(errorstring(Float64(row[i]),Float64(row[i+1]),nsig=1))")
            end
            i < 14 &&  write(io," & ")
        end
        # tex line break
        write(io," \\\\")
        # insert a \hline if the mass changes, insert two of them if beta changes
        if j < nrow
            if data[j,1] != data[j+1,1]
                write(io," \\hline\\hline")
            elseif data[j,2] != data[j+1,2]
                write(io," \\hline")
            end
        end
        # source code line break
        write(io,"\n")
    end

    write(io,"\t\\hline\\hline\n")
    write(io,"\\end{tabular}")
    close(io)
end
write_overview()