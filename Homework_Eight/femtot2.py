from matplotlib import pyplot as plt
import numpy as np


def femtot2():
    data = np.genfromtxt('inputfil.txt',delimiter='')
    n_nodes = int(data[0,0])
    n_cells = int(data[0,1])
    nodes   = data[1:n_nodes+1,0]

    # get x
    x = data[1:n_nodes+1,1]

    
    # get cell data
    cell_start   = n_nodes+1
    cell_end     = cell_start+n_cells
    cells        = data[cell_start:cell_end,0]
    cell_to_node = data[cell_start:cell_end,1:4]
    
    # get node data
    eps_start = cell_end
    eps_end   = eps_start+n_cells
    epsilon   = data[eps_start:eps_end,1]


    # initialize variables
    k0 = 2*np.pi
    n_unknowns = n_nodes
    Z=np.zeros((n_unknowns,n_unknowns),dtype=complex)
    RHS=np.zeros(n_unknowns,dtype=complex)
    


    if False:
        print(f'number of nodes: {n_nodes}')
        print(f'number of cells: {n_cells}')
        # print(f'nodes: {nodes}')
        print(f'cell start: {cell_start}')
        print(f'cell end: {cell_end}')
        print(f'cell_to_node: {cell_to_node}')
        print(f'epsilon start: {eps_start}')
        print(f'epsilon_end: {eps_end}')
        print(f'epsilon: {epsilon}')

    if True:
        print(f'number of nodes: {n_nodes},{np.shape(nodes)}')
        print(f'number of cells: {n_cells},{np.shape(cells)}')
        print(f'number of epsilon: {np.shape(epsilon)}')
        print(f'number of unknowns: {n_unknowns}')



    # fill global matrix one cell at a time using the elemental matrix
    for icell in range(n_cells):
        eleS, eleT = elemat(icell,cell_to_node,x)

        for i in range(3):
            ig = int(cell_to_node[icell,i])
            for j in range(3):
                jg = int(cell_to_node[icell,j])
                Z[ig-1,jg-1]=Z[ig-1,jg-1]+eleS[i,j]-epsilon[icell]*eleT[i,j]*k0**2
    
    # add boundary conditions
    Z[0,0] = Z[0,0] + 1j*k0
    Z[n_nodes-1,n_nodes-1] = Z[n_nodes-1,n_nodes-1] + 1j*k0

    # fill excitation vector
    RHS[0] = 2*1j*k0

    E = np.linalg.inv(Z).dot(RHS)

    plt.figure()
    plt.plot(np.abs(E))
    plt.savefig('E.png')

    gamma = E[0]-1
    mag = np.abs(gamma)
    phs = 180*np.arctan2(np.imag(gamma),np.real(gamma))/np.pi
    print(f'magnitude of reflection: {mag}\n phase of reflection: {phs}\n')


    tau = E[n_unknowns-1]
    mag= np.abs(tau)
    phs = 180*np.arctan2(np.imag(tau),np.real(tau))/np.pi
    print(f'magnitude of transmission: {mag}\n phase of transmission: {phs}')

    




def elemat(icell,cell_to_node, x):
    eleS = np.zeros((3,3))  # derivatives
    eleT = np.zeros((3,3))  # actual functions

    n1 = int(cell_to_node[icell,0])
    n2 = int(cell_to_node[icell,1])
    # n3 = cell_to_node[icell,2]

    delta = x[n2] - x[n1]

    eleS[0,0]=7
    eleS[0,1]=-8
    eleS[0,2]=1
    
    eleS[1,0]=-8
    eleS[1,1]=16
    eleS[1,2]=-8
    
    eleS[2,0]=1
    eleS[2,1]=-8
    eleS[2,2]=7
    

    eleT[0,0]=4
    eleT[0,1]=2
    eleT[0,2]=-1
    
    eleT[1,0]=2
    eleT[1,1]=16
    eleT[1,2]=2
    
    eleT[2,0]=-1
    eleT[2,1]=2
    eleT[2,2]=4

    return eleS/(6*delta), delta*eleT/15
    


if __name__ == '__main__':
    print('in main')
    femtot2()