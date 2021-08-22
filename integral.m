function I = integral(a,b,N)
% Input: %
%           a, b = lower and upper limits of the integral
%           N = number of segments to use in the integration
%           I = approximate value of the integral from a to b of f(x)
%
f = @ (x) sin(x)+x.^ 2; % The function, f(x), we want to integrate

h=(b-a)/N;

x = a:h:b;

x = x';
n = size(x,1);
t = 0;
for i = 2:n
    t = t + h*f((x(i) - x(i-1))/2);
end
I = t;






end