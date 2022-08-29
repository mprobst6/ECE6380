using DelimitedFiles


function mesh1D()
    """

    Generate a 1D finite difference/element mesh1D
    August 29, 2022 M. J. Probst

    """
    n_nodes = 10 # number of n_nodes
    size = 3 # size of the simulation (wavelength)

    delta = size/(n_nodes-1)

    x = range(start=0,step=delta,length=n_nodes)
    e = range(start=1.0,step=0,length=n_nodes)

    open("inputfil.txt","w") do io
        writedlm(io,x)
        writedlm(io,e)
    end
end

mesh1D()