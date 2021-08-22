clear all
format long
%%%%%%%%%%%% Q1 %%%%%%%%%%%%%%%%
k = 3e6;
m1 = 1e4; m2 = 9e3; m3 = 8e3;
A = [2*k/m1, -k/m1, 0; -k/m2, 2*k/m2, -k/m2; 0, -k/m3, k/m3];
eigenval = eigenvals(A);
w = sqrt(eigenval)

%%%%%%%%%%%% Q3 %%%%%%%%%%%%%%%%%%%%

x = 0:24:360;
E = 29e6; I = 720;
y = [0 -0.111 -0.216 -0.309 -0.386 -0.441 -0.473 -0.479 -0.458 -0.412 -0.345 -0.263 -0.174 -0.090 -0.026 0];
delx = 24;
n = size(y,1);
double_deriv = threepointdiff(delx,y')
M_x = E*I.*double_deriv

%%%%%%%%%%%% extra credit question 1 %%%%%%%%%%%%%%%%%%
t = -2100:2100;
c = 4491; a = -2100; b = 2100;
f = (exp(t/c) + exp(-t/c))./2;

trapzoidal_method = trapz(t,f)
F = @(x) (exp(x/c) + exp(-x/c))/2;
Quad_method = quad(F,a,b)

%%%%%%%%%%%%%%%% Q2 Extra credit %%%%%%%%%%%%%%%%%%%%%%%
T = 1e4; E = 29e6; I = 121; w = 1e4;L = 120;
x = 0:30:L;
x = x';
N = size(x,1);
dx = x(2) - x(1);
y = zeros(N,1);
diagT = -(2 + (T*dx^2)/(E*I));
A = eye(N-2,N-2)*diagT;

for i = 1:N-3
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end

C = zeros(N-2,1);
for i = 1:N-2
    C(i) = (w*dx^2)*(L-x(i))*x(i)/(2*E*I);
end
y(2:N-1) = linsolve(A,C);
deflection = y


