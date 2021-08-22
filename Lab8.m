clear all
format long
 %%%%%%%%%%%%% Question 1 B %%%%%%%%%%%%%%%%%
xtol = 1e-5;
ftol = 1e-5;
A = [2,7,5;3,3,4;4,14.001,10];
b = [1;1;1];
x0 = [0;0;0];
maxit = 10;
solA = newtonSys_modified(A, b, x0, xtol, ftol, maxit)

B = [2,7,5;3,3,4;4,14,10];
solB = newtonSys_modified(B, b, x0, xtol, ftol, maxit)

C = [2.1,-6.2,4.3;4.2,-12.5,8.61;-1,-2,3];
solC = newtonSys_modified(C, b, x0, xtol, ftol, maxit)

%%%%%%%%%%%%%%%%%%%%%%%%%% Q 2 B %%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter = 1e5;
xtol = 1e-6;
ftol = 1e-6;
g = @(k) [-110 + k(1)*exp(0.15*k(2)) + 0.15*k(3); -130 + k(1)*exp(0.25*k(2)) + 0.25*k(3); -160 + k(1)*exp(0.35*k(2)) + 0.35*k(3)];

J1 = @(k) [exp(0.15*k(2)), 0.15*exp(0.15*k(2)), 0.15; exp(0.25*k(2)), 0.25*k(1)*exp(0.25*k(2)), 0.25; exp(0.35*k(2)), 0.35*k(1)*exp(0.35*k(2)) , 0.35];
x0 = [10;2;10];
sol_matrix = fsolve(g,x0);
%[v,k,sol] = newtonSys_modified(J1, g, x0, xtol, ftol, max_iter);

k_val = sol_matrix

%%%%%%%%%%% Q 2 C %%%%%%%%%%%%%%%%
F = 100e3; 
r0 = 1;
f = @(r) F/(pi*r^2) - k_val(1)*exp(k_val(2)*r) - k_val(3)*r;
radius = fsolve(f,r0)




