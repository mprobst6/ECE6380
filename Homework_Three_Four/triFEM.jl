using DelimitedFiles
using Plots
using Images

include("conjgrad.jl")

function TriFEM(input_file)

    # read the mesh from the text file
    line = 0

    nnodes::Int64 = 0 # number of nodes
    ncells::Int64 = 0 # number of cells
    ninner::Int64 = 0 
    nouter::Int64 = 0
    nunks::Int64 = 0   

    x = Vector{Float64}() # node positions in x
    y = Vector{Float64}() # node positions in y

    r1 = Vector{Int64}()
    r2 = Vector{Int64}()
    r3 = Vector{Int64}()

    epsilon_r = Vector{Int64}()
    
    f = open(input_file,"r") do f

        # read the parameters
        params::String = readline(f)
        all_params = split(params)

        nnodes = parse(Int64,all_params[1])
        ncells = parse(Int64,all_params[2])
        ninner = parse(Int64,all_params[3])
        nouter = parse(Int64,all_params[4])
        nunks = nnodes - ninner - nouter
        
        for i in 1:nnodes
            line = readline(f)
            values = split(line)
            # print("type of object to be appended: ",parse(Float64,values[2]))
            append!(x,parse(Float64,values[2]))
            append!(y,parse(Float64,values[3]))
        end
        # display(x)
        # print("\n\n")
        # display(y)


        for i in 1:ncells
            line = readline(f)
            values = split(line)
            append!(r1,parse(Int64, values[2]))
            append!(r2,parse(Int64, values[3]))
            append!(r3,parse(Int64, values[4]))
        end
        


        for i in 1:ncells
           line = readline(f)
           values = split(line)
           append!(epsilon_r,parse(Int64,values[2])) 

        end
    end

    # print(size(x))
    pcetond = [r1 r2 r3]
    # Initialize arrays and vectors
    Wtilda = zeros(Float64,nnodes,nnodes)
    W = zeros(Float64,nunks,nunks)
    V = zeros(Float64,nunks)
    elem = zeros(Float64,3,3) # elemental matrix, contribution of one cell to the entries of W
    
    xe = zeros(Float64,3)
    ye = zeros(Float64,3)
    b = zeros(Float64,3)
    c = zeros(Float64,3)

    # loop through the cells, filling the global matrix one cell at a time
    
    for icell in 1:ncells
        n1 = r1[icell]
        n2 = r2[icell]
        n3 = r3[icell]

    # compute 3x3 element matrix
        xe[1] = x[n1]
        xe[2] = x[n2]
        xe[3] = x[n3]

        ye[1] = y[n1]
        ye[2] = y[n2]
        ye[3] = y[n3]

        b[1] = ye[2] - ye[3]
        b[2] = ye[3] - ye[1]
        b[3] = ye[1] - ye[2]

        c[1] = xe[3] - xe[2]
        c[2] = xe[1] - xe[3]
        c[3] = xe[2] - xe[1]

        Area = abs(b[3]*c[1] - b[1]*c[3])*0.5
        # confirmed that the elemental matrices are the same
        for ii = 1:3
            for jj = 1:3
                elem[ii,jj] = (b[ii]*b[jj] + c[ii]*c[jj])*epsilon_r[icell]/Area/4
            end
        end

       
        for ii in 1:3
            ig = pcetond[icell,ii] # ig is the global node for ii
            for jj in 1:3
                jg = pcetond[icell,jj] # jg is the global node for jj
                Wtilda[ig,jg] = Wtilda[ig, jg] + elem[ii,jj]
                if (ig <= nunks)
                    if (jg <= nunks)
                        W[ig,jg] = W[ig,jg] + elem[ii,jj]
                    elseif (jg <= (nunks + nouter))
                        V[ig] = V[ig] - elem[ii,jj]
                    end
                end
            end
        end
    end


    CG_Pot = ones(length(V))
    CG_Pot = conjgrad(W,V,CG_Pot)

    Pot = W\V


    Pot = [Pot; ones(nouter); zeros(ninner)]
    CG_Pot = [CG_Pot; ones(nouter); zeros(ninner)]
     
    Cap = (Pot')*Wtilda*Pot

    CG_Cap = (CG_Pot')*Wtilda*Pot

    str = string("Capacitance/epsilon0 pul = ",Cap,"\n")
    CG_str = string("Capacitance/epsilon w/conjgrad = ",CG_Cap,"\n")
    ground_truth = 2*pi/(log(2.5))
    
    # print(ncells,"\n")
    print(str)
    print(CG_str)
    print("Error of conjugate gradient: ",(1.0-(CG_Cap/ground_truth))*100,"%\n")

    open("potfil_jl.txt","w") do io
        write(io,str)  
        str = "node   Potential\n"
        write(io,str)

        for ii in 1:nunks
            str = string(ii,"  ",real(Pot[ii]),"\n")
            write(io,str)
        end
    end

    return
end


input_file = "cylfil.txt"
TriFEM(input_file)
