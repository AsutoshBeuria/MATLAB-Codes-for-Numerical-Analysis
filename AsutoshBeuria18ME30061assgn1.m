
close all
syms x;
f1 = @(x) 3*x + 4; g1 = 0;
f2 = @(x) sin(2*pi*x); g2 = @(x) -(4*pi^2)*sin(2*pi*x);
f3 = @(x) x.^3 + x; g3 = @(x) 6*x;

dx = [0.1; 0.05; 0.01; 0.001];

h1(x) = secondorder(f1,dx,x)
h2(x) = secondorder(f2,dx,x)
h3(x) = secondorder(f3,dx,x)

p1(x) = thirdorder(f1,dx,x)
p2(x) = thirdorder(f2,dx,x)
p3(x) = thirdorder(f3,dx,x)

r1(x) = fourthorder(f1,dx,x)
r2(x) = fourthorder(f2,dx,x)
r3(x) = fourthorder(f3,dx,x)


function y = secondorder(f,dx,x)
y = (f(x+dx) - 2*f(x) + f(x-dx))./(dx.^2);
y = simplify(y);
y = vpa(y,5);
end

function y = thirdorder(f,dx,x)
a = 35/12; b = -26/3;
c = 19/2; d = -14/3; 
e = 11/12;
y = (a*f(x) + b*f(x+dx) + c*f(x+ 2*dx) + d*f(x+ 3*dx) + e*f(x + 4*dx))./(dx.^2);
y = simplify(y);
y = vpa(y,5);
end

function y = fourthorder(f,dx,x)
a = -1/12; b = 16/12;
c = -30/12; d = 16/12;
e = -1/12;

y = (a*f(x + 2*dx) + b*f(x+dx) + c*f(x) + d*f(x - dx) + e*f(x - 2*dx))./(dx.^2);
y = simplify(y);
y = vpa(y,5);
end


