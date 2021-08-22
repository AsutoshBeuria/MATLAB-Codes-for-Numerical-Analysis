clear all
format long
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





