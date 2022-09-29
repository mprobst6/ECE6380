function power()
    A = [3 1; 1 2]
    e = [1; 1]
    e_new = Vector{Float64}()

    for i in 1:7
        e_new = A*e
        
        lambda::Float64 = maximum(e_new)
        e_new = e_new/lambda
        print("\nlambda: ",lambda)
        print("\nnormalized e_new: ",e_new,"\n")
        e = e_new
        if mod(i,5) == 0
            print("approximate eigenvalue: ",lambda,"\n")
        end
    end
    return
end

power()
