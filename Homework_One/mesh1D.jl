using DelimitedFiles


function mesh1D()
    """
    Generate a 1D finite difference/element mesh1D
    August 29, 2022 M. J. Probst

    """
    input_file = string("inputfil",ARGS[1],".txt")
    n_nodes = parse(Int16, ARGS[1]) # number of n_nodes
    domain = 1.5 # size of the simulation (wavelength)

    delta = domain/(n_nodes-1)

    x = range(start=0,step=delta,length=n_nodes)
    e = range(start=1.0,step=0,length=n_nodes)

    open(input_file,"w") do io
        writedlm(io,[x e],',')
    end
end

mesh1D()