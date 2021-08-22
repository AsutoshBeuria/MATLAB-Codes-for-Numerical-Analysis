clear all
format long
F = @(x) [x(1) + x(2) - 2; x(1)*x(3) + x(4)*x(2); x(1)*x(3)^2 + x(2)*x(4)^2 - (2/3); x(1)*x(3)^3 - x(2)*x(4)^3];
x0 = [1;1;0.1;0.1];
[x,fval] = fsolve(F,x0)

x0 = [1;-1;-0.1;0.1];
[x,fval] = fsolve(F,x0)

x0 = [-2;1;0.1;1];
[x,fval] = fsolve(F,x0)

x0 = [1;1;0;0.1];
[x,fval] = fsolve(F,x0)

x0 = [0;0;1;3];
[x,fval] = fsolve(F,x0)

%%%%%%% Q2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun = @(x) [(x(1) + 3)*(x(2)^3 - 7) + 18;sin(x(2)*exp(x(1)) - 1)];

J = @(x) [x(2).^3 - 7.0, 3*x(2).^2*(x(1) + 3) ; x(2)*exp(x(1))*cos(x(2)*exp(x(1)) - 1), exp(x(1))*cos(x(2)*exp(x(1)) - 1)];

x0 = [-0.5;1.4]; 
max_iter = 1e5;
xtol = 1e-6;
ftol = 1e-6;

[v,k,sol] = newtonSys_xx(J, fun, x0, xtol, ftol, max_iter);

solution = sol
val_of_xtillconvergence = v
x1_val = v(:,1);
x2_val = v(:,2);
trueerr_x1 = abs(x1_val - 0);
trueerr_x2 = abs(x2_val - 1);

no_of_iter = [1:k];

data_table1 = table(no_of_iter', x1_val, trueerr_x1)
data_table2 = table(no_of_iter', x2_val, trueerr_x2)





