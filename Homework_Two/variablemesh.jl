using DelimitedFiles


function variable_mesh(array,file)
    """
    Generate a 1D finite difference/element mesh1D with variable spacings
    August 29, 2022 M. J. Probst

    """
    folder = "inputfiles"
    mkpath(folder)
    input_file = string("inputfiles/",file,".txt")

    open(input_file,"w") do io
        writedlm(io,[array[:,1] array[:,2]],',')
    end
end

# code to generate the input mesh stored in inputfil_1.txt
if true # true
    n_nodes::Int8 = 25
    mesh_idxs = range(start=1,stop=n_nodes,length=n_nodes)
    eps_idxs = range(start=1,stop=24,length=n_nodes-1)
    mesh_vals = vec([0 0.05 0.1 0.15 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5 0.525 0.55 0.575 0.6 0.65 0.7 0.75 0.8])
    eps_vals = [range(start=1,stop=1,length=4); range(start=5, stop=5, length=16); range(start=1,stop=1,length=4)]

    input_array = [n_nodes 0; mesh_idxs mesh_vals; eps_idxs eps_vals]

    variable_mesh(input_array,"inputfil_1")
end


# code to generate the input mesh stored in inputfil_2.txt
if true # true
    n_nodes::Int8 = 29
    mesh_idxs = range(start=1, stop=n_nodes,length=n_nodes)
    eps_idxs = range(start=1,stop=n_nodes-1,length=n_nodes-1)
    mesh_vals = vec([0 0.05 0.1 0.15 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.52 0.54 0.56 0.58 0.6 0.65 0.7 0.75 0.8])
    eps_vals = [range(start=1,stop=1,length=4); range(start=5,stop=5,length=20); range(start=1,stop=1,length=4)]


    # try without a variable mesh...if I get closer to their thing or magnitude
    # mesh_vals = range(start=0,stop=0.8,length=n_nodes)
    # eps_vals = [range(start=1,stop=1,length=20); range(start=5,stop=5,length=40); range(start=1,stop=1,length=20)]

    input_array = [n_nodes 0; mesh_idxs mesh_vals; eps_idxs eps_vals]

    variable_mesh(input_array,"inputfil_2")
end
