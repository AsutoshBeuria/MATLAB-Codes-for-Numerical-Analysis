tspan = [0 10]; y0 = 0;
a = 1; b = 1; A = 3; g = 10; rho = 0.1;
[t,y] = ode45(@(t,y) (a/A - (a*rho*g*y/(b*A))), tspan, y0);          % Actually in the previus code that i had sent you , i wrote the expression as (2*a/A - (a*rho*g*y/(b*A)))
plot(t, y, '-o')                                                       %So I think thats why I got wrong plot.

xlabel('time(s) ') ; ylabel('height(m)') ; title('h(t) graph') 

