include("conjgrad.jl")

size = 20
# A = rand(size,size)
# b = rand(size,1)
# x = rand(size,1)
A = [1 1 1 1; 2 4 6 8; 1 4 9 16; 4 3 2 1]
b = [1; 5; 3; 6]
x = [1; 0; 0; 0]

# print("A: ",A,"\n\n")
# print("b: ",b,"\n\n")
x = conjgrad(A,b,x)
# print("final x: ",x,"\n\n")
# print("final b: ",b,"\n\n")
# print("A*x: ",A*x,"\n\n")