using DelimitedFiles
using Plots

# modify the code to take in the numebr of nodes at teh start
# make sure I have 25 points and only 24 refractive indices

function setup_sim(input_file)
    """
    Read the input file and setup parameters necessary for simulation
    """
    # read data
    f = open(input_file,"r")
    data = readdlm(input_file, ',', Float64)

    n_nodes::Int8 = data[1,1]
    nodes=data[2:n_nodes+1,2]
    epsilon = data[n_nodes+2:end,2]



    # initialize variables
    k0::Float32 = 2*pi
    Z = zeros(Complex{Float32},n_nodes,n_nodes)
    R = zeros(Complex{Float32},n_nodes)

    sim_stuff = Dict(
        "nodes"   => nodes,
        "epsilon" => epsilon,
        "Z"       => Z,
        "R"       => R
    )

    params = Dict(
        "n_nodes" => n_nodes,
        "k0"      => k0
    )
    return params, sim_stuff
end


function scattering(params, sim_stuff)
    """
    Fill the global matrix for scattering fields with variable mesh sizes
    """
    k0 = params["k0"]
    n_nodes = params["n_nodes"]

    nodes = sim_stuff["nodes"]
    epsilon = sim_stuff["epsilon"]
    Z = sim_stuff["Z"]
    R = sim_stuff["R"]

    # print(size(epsilon))
    # print(size(Z))
    # print(size(R))
    # print(size(nodes))

    for row in 1:n_nodes
        if row == 1
            # define difference sizes
            delta_R = nodes[row+1]-nodes[row]
            print("\nrow",row,"   delta_R: ",round(delta_R; digits=3))

            # fill global matrix
            Z[row,row] = 1/delta_R - (k0^2)*delta_R*epsilon[row]/3 + im*k0
            Z[row,row+1] = -1/delta_R - (k0^2)*delta_R*epsilon[row]/6

            # fill excitation vector (right hand side)
            R[row] = im*2*k0
        elseif row == n_nodes
            # define difference sizes
            delta_L = nodes[row] - nodes[row-1]
            print("\nrow",row,"   delta_L: ",round(delta_L;digits=3))

            # fill global matrix
            Z[row,row-1] = -1/delta_L - (k0^2)*delta_L*epsilon[row-1]/6
            Z[row,row] = 1/delta_L - (k0^2)*delta_L*epsilon[row-1]/3 + im*k0
        else
            # define difference sizes
            delta_R = nodes[row+1] - nodes[row]
            delta_L = nodes[row] - nodes[row-1]
            print("\nrow",row,"   delta_L: ", round(delta_L;digits=3), "  and   delta_R: ",round(delta_R;digits=3))


            # fill global matrix
            Z[row,row-1] = -1/delta_L - (k0^2)*delta_L*epsilon[row-1]/6
            Z[row,row] = 1/delta_L + 1/delta_R - (k0^2)*( (epsilon[row-1]*delta_L) + (epsilon[row]*delta_R) )/3
            Z[row,row+1] = -1/delta_R - (k0^2)*epsilon[row]*delta_R/6
        end
    end

    # solve the system of equations for the field data
    E = Z \ R 
    return E
end

function save_data(field, sim_stuff, coding_output)
    """
    Output the transmission and reflection data to an output file and plot certain quantities
    """

    mkpath("outputfiles") # used to store all output files

    file = string("output")
    nodes = sim_stuff["nodes"]
    epsilon = sim_stuff["epsilon"]

    # Extract the real, imaginary, and absolute values
    R = Array{Float32,1}(undef, length(nodes))
    I = Array{Float32,1}(undef, length(nodes))
    A = Array{Float32,1}(undef, length(nodes))

    # fill the vectors
    for i in 1:length(nodes)
        R[i] = real(field[i])
        I[i] = imag(field[i])
        A[i] = abs(field[i])
    end
    P = [R, I, A]
    plot(nodes, P, label=["Real" "Imag" "Abs"],lw=3)
    savefig(string(file,".png"))

    plot(nodes)
    savefig("nodes.png")

    plot(nodes[1:end-1],epsilon)
    savefig("epsilon.png")

    gamma = field[1] - 1
    mag = abs(gamma)
    phs = 180*atan(imag(gamma), real(gamma))/pi

    
    ref_mag = string("Magntiude: ",mag,"\n")
    ref_phase = string("Phase: ",phs,"\n\n")

    tau = field[end]
    mag = abs(tau)
    phs = 180*atan(imag(tau),real(tau))/pi
    tran_mag = string("Magnitude: ",mag,"\n")
    tran_phase = string("Phase: ",phs,"\n")

    open(coding_output,"w") do io
        write(io,"Reflection\n")
        write(io, ref_mag)
        write(io, ref_phase)
        write(io, "Transmission\n")
        write(io, tran_mag)
        write(io, tran_phase)
    end
    return gamma, tau
end




##########################################
### Programming Assignment Problem One ###
##########################################

# input file 1
if true
    coding_input = "inputfiles/inputfil_1.txt"
    coding_output = "outputfiles/outputfil_1.txt"
    
    coding_params, coding_data = setup_sim(coding_input)
    E = scattering(coding_params, coding_data)
    save_data(E, coding_data,coding_output)
end

# input file 2
if true
    coding_input = "inputfiles/inputfil_2.txt"
    coding_output = "outputfiles/outputfil_2.txt"
    
    coding_params, coding_data = setup_sim(coding_input)
    E = scattering(coding_params, coding_data)
    save_data(E, coding_data,coding_output)
end



######################################
### Regular Assignment Problem One ###
######################################
if false
    gamma = Array{Complex{Float64},1}(undef,6)
    tau = Array{Complex{Float64},1}(undef,6)
    
    for i in 2:6
        interp_input = string("inputfiles/slab_",2^i,".txt")
        interp_output = string("outputfiles/slab_",2^i,".txt")
        params, sim_stuff = setup_sim(interp_input)
        E = scattering(params, sim_stuff)
        

        gamma[i-1],tau[i-1] = save_data(E, sim_stuff, interp_output)
    end

    gamma_exact = Array{Complex{Float64},1}(undef,5)
    tau_exact = Array{Complex{Float64},1}(undef,5)
    ref_mag = Array{Float64,1}(undef,5)
    ref_phs = Array{Float64,1}(undef,5)
    tran_mag = Array{Float64,1}(undef,5)
    tran_phs = Array{Float64,1}(undef,5)

    for i in 1:5
        # calculate the "exact" values
        gamma_exact[i] = (4*gamma[i+1] - gamma[i])/3
        tau_exact[i] = (4*tau[i+1] - tau[i])/3

        # calculate the magnitude and phase of the "exact" values
        ref_mag[i] = abs(gamma_exact[i])
        ref_phs[i] = 180*atan(imag(gamma_exact[i]), real(gamma_exact[i]))/pi
        tran_mag[i] = abs(tau_exact[i])
        tran_phs[i] = 180*atan(imag(tau_exact[i]), real(tau_exact[i]))/pi
    end
    
    final_output = "outputfiles/interpolate.txt"
    open(final_output,"w") do io
        write(io,"Reflection Magnitudes\n")
        for i in 1:4
            write(io,string(2^(i+1),"/",2^(2+i)," cells: ",ref_mag[i],"\n"))
        end
        write(io,"\nReflection Phases\n")
        for i in 1:4
            write(io,string(2^(i+1),"/",2^(2+i)," cells: ",ref_phs[i],"\n"))
        end

        write(io,"\n\n\nTransmission Magnitudes\n")
        for i in 1:4
            write(io,string(2^(i+1),"/",2^(2+i)," cells: ",tran_mag[i],"\n"))
        end
        write(io,"\nTransmission Phases\n")
        for i in 1:4
            write(io,string(2^(i+1),"/",2^(2+i)," cells: ",tran_phs[i],"\n"))
        end
    end 
end
