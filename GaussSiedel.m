clear all
close all
format long 
A = [4,0,1,0,1;2,5,-1,1,0;1,0,3,-1,0;0,1,0,4,-2;1,0,-1,0,5];
b = [32;19;14;-2;41];
x0 = [1;1;1;1;1];

n_max = 1e5;
tol = 1e-12;

[X_iter,sol_GS,resid,relerr,iterno_GS] = GaussMethod(A,b,n_max,x0,tol)

true_sol = A\b;
no_iterationsJacobi = Jacobi(A,b,n_max,x0,tol);
x = linspace(1,iterno_GS,size(relerr,1));

plot(x, relerr,'-o','Linewidth',2);
title('Convergence plot-Norm P2 error'); xlabel('No. of iterations'); ylabel('Error');
figure();

plot(x, resid,'-o','Linewidth',2);
title('Convergence plot- Residual'); xlabel('No. of iterations'); ylabel('Error');




