clear all
A = [8,-1,0,4;2,0,1,1;8,15,-2,0;1,4,6,-13]
b = [1;-1;10;4]
n = size(A,1);

[L1,U1,P1] = LUPDecompos(A,b)

[L2,U2] = Croutmethod(A)

[L3,U3,P3] = lu(A,'matrix')

y = zeros(n,1); x = zeros(n,1);  % LU = A, Ax = b; This becomes LUx = b;, Now we have to solve for y = Ux And Ly = b
y(1) = b(1)/L2(1,1); 
% Forward Substitution                 
for i = 2:n
    y(i) = (b(i) - L2(i,1:i-1)*y(1:i-1))./L2(i,i);
end
% Backward Substitution
x(n) = y(n);
for i = n-1:-1:1      
    x(i) = (y(i) - U2(i,i+1:n)*x(i+1:n));
end
sol_LU = x

n_max = 1e5;
x0 = ones(n,1);
tol = 0.01;
[sol_GS, P_GS, no_iter] = GaussModified(A,b,n_max,x0,tol)






