clear all;
S = [1,0,-0.5;0,1,-0.5;-0.5,-0.5,1];

E = 200e9; 
A = 8e-03; 
L = 2; P = 20e3;
B = (E*A/L)*S;
c = [0;0;-P];

mNAM_b = [c B(:,2:3)]
delta_0 = det(B)
delta_b = det(mNAM_b) 
VB = delta_b/delta_0
% Create a modified NAM for node c by inserting I into the second column.
mNAM_c = [B(:,1) c B(:,3)]
delta_c = det(mNAM_c)
VC = delta_c/delta_0
% Create a modified NAM for node e by inserting I into the third column.
mNAM_e = [B(:,1) B(:,2) c]
delta_e = det(mNAM_e)
VE = delta_e/delta_0

Gauss(B,c)