function FEM

    global xy cell_to_node cell_to_edge epsilon edge_to_cell edge_to_node;

    nnodes = dlmread('cylfil.txt','',[0,0,0,0]);
    ncells = dlmread('cylfil.txt','',[0,1,0,1]);
    nedges = dlmread('cylfil.txt','',[0,2,0,2]);
    ninner = dlmread('cylfil.txt','',[0,3,0,3]);

    xy = dlmread('cylfil.txt','',[1,1,nnodes,2]);

    nstart = nnodes+1;
    nend = nstart + ncells -1;
    cell_to_node = dlmread('cylfil.txt','',[nstart,1,nend,3]);

    nstart = nend + 1;
    nend = nstart + ncells - 1;
    cell_to_edge = dlmread('cylfil.txt','',[nstart,1,nend,3]);

    W = zeros(nedges);
    Y = zeros(nedges);
   

    for icell=1:ncells
        [U,V] = elemat(icell); % need to write this
        
        n1 = cell_to_node(icell,1); n2 = cell_to_node(icell,2); n3 = cell_to_node(icell,3);
        
        % Orient the basis vectors so that the point inward
        for kk=1:3
            iedge=0;
            if ((kk==1) && (n2 > n3))
                iedge = 1;
            end
            if ((kk==2) && (n3 > n1))
                iedge = 1;
            end
            if ((kk==3) && (n1 > n2))
                iedge = 1;
            end
            if (iedge == 1)
                for jj = 1:3
                    V(kk,jj) = -V(kk,jj);
                    V(jj,kk) = -V(jj,kk);
                    U(kk,jj) = -U(kk,jj);
                    U(jj,kk) = -U(jj,kk);
                end
            end
        end
  
        % Fill the global matrices
        for jj = 1:3
            jp = cell_to_edge(icell,jj);
            for kk = 1:3
                kp = cell_to_edge(icell,kk);
                W(jp,kp) = W(jp,kp) + U(jj,kk);
                Y(jp,kp) = Y(jp,kp) + V(jj,kk);
            end
        end
    end
    
    fid = fopen('eigfil.txt','wt');

    E = sort(eig(W,Y));
    str = 'TM resonant waveneumbers: ';
    fprintf(fid,'%s \n',str);
    for ii = 1:nedges
        reaE = real(sqrt(E(ii)));
        fprintf(fid,'%6d %15.14g\n',ii,reaE);
    end

end
% -------------------------------------------------------------------------
function [U,V] = elemat(icell)
    global xy cell_to_node;
    U(3,3) = 0;
    V(3,3) = 0;

    % Need the node locations to find the side edges etc
    n1 = cell_to_node(icell,1); n2 = cell_to_node(icell,2); n3 = cell_to_node(icell,3);

    x(1) = xy(n1,1); y(1) = xy(n1,2);
    x(2) = xy(n2,1); y(2) = xy(n2,2);
    x(3) = xy(n3,1); y(3) = xy(n3,2);

    b(1) = y(2)-y(3); b(2) = y(3)-y(1); b(3) = y(1)-y(2);
    c(1) = x(3)-x(2); c(2) = x(1)-x(3); c(3) = x(2)-x(1);

    Area = abs(b(3)*c(1) - b(1)*c(3))*0.5;

    % Pythagorean Theorem to find the length of the sides
    w(1) = sqrt((x(2)-x(3))^2 + (y(3)-y(2))^2);
    w(2) = sqrt((x(3)-x(1))^2 + (y(3)-y(1))^2);
    w(3) = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);

    for m = 1:3
        for n = 1:3
            U(m,n) = w(m)*w(n)*(b(mod(m,3)+1)*c(mod(m-2,3)+1)-b(mod(m-2,3)+1)*c(mod(m,3)+1))*(b(mod(n,3)+1)*c(mod(n-2,3)+1)-b(mod(n-2,3)+1)*c(mod(n,3)+1));
        end
    end

    U = U/(4*Area^3);

    for m = 1:3
        for n = 1:3
            for i = 1:2
                for j = 1:2
                    if (i==j)
                        alpha = 1;
                    else
                        alpha = -1;
                    end
                    if (m+i == n+j)
                        gamma = 1/12;
                    else
                        gamma = 1/24;
                    end
                    V(m,n) = V(m,n) + alpha*gamma*(b(mod(m+2-i,3)+1)*b(mod(n+2-j,3)+1)+c(mod(m+2-i,3)+1)*c(mod(n+2-j,3)+1));
                end
            end
            V(m,n) = V(m,n)*w(m)*w(n);
        end
    end
    V = V/(2*Area);
    % disp(U);
    % disp(V);
end

0.94079249422024
0.98804419
0.995002


average edge lengths: 1.8340 0.7139 0.4438