using DelimitedFiles
function slab(cells)
    """
    Generate a 1D finite element mesh with constant spacing such that the 
    number of cells in the mesh is a parameters
    """
    folder = "inputfiles/"
    mkpath(folder)
    input_file = string(folder,"slab_",cells,".txt")

    n_nodes::Int8 = cells + 1

    # locations of the nodes
    mesh_idx = range(start=1,stop=n_nodes,length=n_nodes)
    mesh = range(start=0,stop=0.4,length=n_nodes)

    
    # epsilon is a slab with epsilon=4 and thickness 0.2lambda, there is 0.1 lambda of air on either side
    eps_idx = range(start=1,stop=Int8(cells),length=Int8(cells))
    eps = [range(start=1,stop=1,length=Int8(cells/4)); range(start=4,stop=4,length=Int8(cells/2)); range(start=1,stop=1,length=Int8(cells/4))]

    array = [n_nodes 0; mesh_idx mesh; eps_idx eps]

    open(input_file,"w") do io
        writedlm(io,[array[:,1] array[:,2]],',')
    end

    return
end

powers = range(start=2, step=1, stop=6)
for power in powers
    n_cells::Int16 = 2^power
    slab(n_cells)
end