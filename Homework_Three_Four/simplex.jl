using LinearAlgebra

function simplex_from_cartesian(x,y)
    """
    Function to calculate simplex coordinates from cartesian coordinates

    """
    # load the vertices
    global vertices
    a1 = vertices[2,1]*vertices[3,2] - vertices[3,1]*vertices[2,2]
    b1 = vertices[2,2]-vertices[3,2]
    c1 = vertices[3,1]-vertices[2,1]

    a2 = vertices[3,1]*vertices[1,2] - vertices[1,1]*vertices[3,2]
    b2 = vertices[3,2]-vertices[3,2]
    c2 = vertices[1,1]-vertices[3,1]

    a3 = vertices[1,1]*vertices[2,2] - vertices[2,1]*vertices[1,2]
    b3 = vertices[1,2]-vertices[2,2]
    c3 = vertices[2,1]-vertices[1,1]

    # find the area
    M = [1 1 1; vertices[:,1]'; vertices[:,2]';] # print this to ensure correctness
    A2 = det(M)

    # calculate the simplex coordinates
    L1 = (a1 + b1*x + c1*y)/A2
    L2 = (a2 + b2*x + c2*y)/A2
    L3 = (a3 + b3*x + c3*y)/A2
    return L1, L2, L3
end


function cartesian_from_simplex(L1,L2,L3)
    """
    Function to calculate cartesian coordinates from simplex coordinates
    """
    global vertices
    x = L1*vertices[1,1] + L2*vertices[2,1] + L3*vertices[3,1]
    y = L1*vertices[1,2] + L2*vertices[2,2] + L3*vertices[3,2]
    return x, y
end


### MAIN ###
# PART A
L1::Float64 = 0.33333333333
L2::Float64 = 0.33333333333
L3::Float64 = 0.33333333333
vertices = [0 0; 0.5 1; 1 0]
x,y = cartesian_from_simplex(L1, L2, L3)
print("x: ",x,"\ny: ",y,"\n")

# PART B
I1::Float64 = 0.75
I2::Float64 = 0.5
O1, O2, O3 = simplex_from_cartesian(I1,I2)
print("\nL1: ",O1,"\nL2: ",O2,"\nL3: ",O3,"\n")