clear all;
S = [1,0,-0.5;0,1,-0.5;-0.5,-0.5,1];

E = 200e9; 
A = 8e-03; 
L = 2; P = 20e3;
K = (E*A/L)*S;
F = [0;0;-P];

[L,U] = LUdecompos(K,F)

