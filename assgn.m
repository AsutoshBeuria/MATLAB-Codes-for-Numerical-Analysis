format long
syms x
f = 1e-03*(-0.900171*x^3 +  4.80091*x^2 - 11.3799);
g(x) =  vpa(diff(f,x),4);