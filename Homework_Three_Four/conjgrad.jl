function conjgrad(A,b,x)
    tolerance::Float64 = 1e-5
    denom::Float64 = only(b' * b)

    r = A*x - b

    p = -r

    rsold::Float64 = only(r' * r) # 34

    rnrm::Float64 = 1

    for i in 1:length(b)
        Ap = A*p

        alpha::Float64 = rsold / only(p' * Ap)

        x = x + alpha*p
 
        r = r + alpha*Ap

        rsnew::Float64 = only(r' * r )

        rnrm = sqrt(rsnew/denom)


        if rnrm < tolerance
            print("Took ",i," iterations\n\n")
            break
        end
        p = -r + (rsnew / rsold)*p

        rsold = rsnew

    end
    return x
end