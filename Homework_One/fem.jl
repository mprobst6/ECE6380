using DelimitedFiles
using Printf
using Plots


# TODO: put the params (nodes, epsilon, delta) in a dictionary and pass that out

function setup_sim(input_file)
    """
    Read the input file and setup parameters necessary for simulation
    """
    # read data
    f = open(input_file,"r")
    data = readdlm(input_file, ',', Float64)

    # initialize variables
    k0::Float32 = 2*pi    # k vector (normalized to wavelength)
    x = data[:,1]         # grid points
    nodes = length(x)     # number of nodes
    epsilon = data[:,2]   # refractive index sampled at each node
    delta = x[2] - x[1]   # spacing between nodes
    Z = zeros(Complex{Float32},nodes,nodes) # initialize system of equations matrix
    R = zeros(Complex{Float32},nodes)       # initialize initial conditions vector

    params = Dict(
        "nodes" => nodes,
        "delta" => delta,
        "k0"    => k0
    )

    return x, epsilon, Z, R, params
end

function solve(Z, R)
    """
    Solve the system of equations whose solution is the electric field
    """
    E = Z \ R  # matrix solve
    return E
end


function dirichlet(x,epsilon, Z, R, params)
    """
    Construct the Z matrix subject to Dirichlet boundary conditions
    """
    delta = params["delta"]
    k0 = params["k0"]
    nodes = params["nodes"]
    
    # specify Dirichlet boundary conditions
    Ka::Complex{Float32} = 1
    Kb::Complex{Float32} = exp(-im*k0*x[end])
    R[1] = Ka
    R[end] = Kb

    # construct the RHS and Z matrices
    global Z[1,1] = 1
    global Z[nodes,nodes] = 1

    for irow = 2:nodes-1
        global Z[irow,irow-1] = 1/delta^2 + k0^2 * epsilon[irow-1]/6
        global Z[irow,irow] = -2/delta^2 + k0^2 * (epsilon[irow-1] + epsilon[irow+1]) / 3
        global Z[irow,irow+1] = 1/delta^2 + k0^2 * epsilon[irow+1]/6
    end
    return Z, R
end


function neumann(x, epsilon, Z, R, params)
    """
    Construct the Z matrix subject to a Neumann boundary condition and a Dirichlet boundary condition
    """

    delta = params["delta"]
    k0 = params["k0"]
    nodes = params["nodes"]

    # speficy the known boundary conditions
    Ka::Complex{Float32} = 1
    R[1] = Ka

    # construct the Z matrix
    Z[1,1] = 1

    for irow = 2:nodes-1
        Z[irow,irow-1] = 1/delta^2 + k0^2 * epsilon[irow-1]/6
        Z[irow,irow] = -2/delta^2 + k0^2 * (epsilon[irow-1] + epsilon[irow+1]) / 3
        Z[irow,irow+1] = 1/delta^2 + k0^2 * epsilon[irow+1]/6
    end

    Z[nodes,nodes-1] = 1/delta^2 + k0^2 * (epsilon[nodes-1])/6
    Z[nodes,nodes] = -1/delta^2 + k0^2 * (epsilon[nodes-1])/3 - im*k0/delta
    return Z, R
end


function save_data(x,E,Z)
    """
    Save the field data (real, imaginary, and magnitude) and plot versus wavelength
    """

    file = string("output_",ARGS[1])

    R = Array{Float32,1}(undef,length(x))
    I = Array{Float32,1}(undef,length(x))
    A = Array{Float32,1}(undef,length(x))

    for i in 1:length(R)
        R[i] = real(E[i])
        I[i] = imag(E[i])
        A[i] = abs(E[i])
    end

    P = [R, I, A]
    p1 = plot(x,E, title="Electric Complex",label="",lw=3)
    p2 = plot(x,P,label=["Real" "Imag" "Abs"],lw=3)
    plot(p1,p2,layout=(2,1),legend=true)
    savefig(string(file,".png"))

    open(string(file,".txt"),"w") do io
        writedlm(io,Z[1:5,1:5])
    end
end



input_file = string("inputfil",ARGS[1],".txt")    # input file (used to iterate over different resolutions)

x, epsilon, Z, R, params = setup_sim(input_file)  # setup the data (boundary condition independent)
Z, R = neumann(x,epsilon, Z, R, params)           # construct Z and R from the given boundary conditions
# Z, R = dirichelet(x, epsilon, Z, R, params)
E = solve(Z, R)                                   # solve the system of equations for the field data
save_data(x, E, Z)                                # save the field data