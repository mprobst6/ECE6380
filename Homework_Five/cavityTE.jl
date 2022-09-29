using LinearAlgebra

function cavityTE(input_file)

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
    W = zeros(Float64,nnodes,nnodes)
    Y = zeros(Float64,nnodes,nnodes)
    eleS = zeros(Float64,3,3) # elemental matrix, contribution of one cell to the entries of W
    eleT = zeros(Float64,3,3)

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
                eleS[ii,jj] = (b[ii]*b[jj] + c[ii]*c[jj])*epsilon_r[icell]/Area/4
                if ii == jj
                    eleT[ii,jj] = Area/6
                else
                    eleT[ii,jj] = Area/12
                end
            end
        end
    # add contributions from cell 'cell' to global matrix'

        for jj in 1:3
            jp = pcetond[icell,jj]
            for kk in 1:3
                kp = pcetond[icell,kk]
                W[jp,kp]=W[jp,kp]+eleS[jj,kk]
                Y[jp,kp]=Y[jp,kp]+eleT[jj,kk]
            end
        end
    end
    
    vals, vects = eigen(W,Y)

    str = "TE resonant wavenumbers: \n"
    print(str)

    for ii in 2:7
        realE = sqrt(vals[ii])
        # realE = real(sqrt(vals[ii]))
        # imagE = imag(sqrt(vals[ii]))
        print("real E: ",realE,"\n") # "  imaginary E: ",imagE,"\n")
    end
    print("\n")
    return
end

input_file = "cylfil.txt"
cavityTE(input_file)