
function y = conjgrad(A, b, x)
    tolerance = 1.0e-5;
    denom = b' * b;
    % disp(denom)
    r = A * x - b;
    p = -r;
    rsold = r' * r;   % 34

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r + alpha * Ap;
        rsnew = r' * r;
        rnrm = sqrt(rsnew/denom)
        y = x;
        if rnrm < tolerance
              break;
        end
        p = -r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end