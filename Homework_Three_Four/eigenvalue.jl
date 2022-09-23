using LinearAlgebra

A = [1 10000;0 2]
B = A'
AtA = (B) * (A)
C = eigvals(AtA)
print(B,"\n")
e1 = maximum(C)
e2 = minimum(C)

print("condition number: ",sqrt(e1/e2),"\n")
print("eigenvalues",C,"\n")
