format long
syms x
f = @(x) 1e-03*(-0.900171*x^3 +  4.80091*x^2 - 11.3799);
g(x) =  vpa(diff(f,x),4);

i = 1;
while 1
    i = i - f(i)/g(i);
    y = i;
    if abs(g(y)) < 1e-5
        sol = y
        break
    end  
end