function jacob = J()
syms f1 f2 x1 x2

f1(x1,x2) = (x1 + 3)*(x2^3 - 7) + 18;
f2(x1,x2) = sin(x2*exp(x1) - 1);
fun = 

J(x1,x2) = [diff(f1,x1), diff(f2,x1); diff(f2,x1), diff(f2,x2)];
jacob = matlabFunction(J(x1,x2));