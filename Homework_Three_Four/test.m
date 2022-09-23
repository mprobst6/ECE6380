A = [1 1 1 1; 2 4 6 8 ; 1 4 9 16; 4 3 2 1];
b = [1; 5; 3; 6];
x = [1; 0; 0; 0];
y = conjgrad(A,b,x);
% disp(y)
% disp(A*y)