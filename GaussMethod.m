function [P,sol,resid,relerr,iterno] = GaussMethod(A,b,n_max,x0,tol) 
% Iteration k = 1 ;
x = x0; %w = 1.045; % min w value obtained from numerical experiment
V = [x];
resid = [];
relerr = [];
n = size(A,1);
for iter = 1:n_max 
    x_old = x;
    for i = 1:n
        sigma = 0;
        for j = 1:i-1
            sigma=sigma+A(i,j)*x(j);
        end
        for j=i+1:n
                sigma=sigma+A(i,j)*x_old(j);
        end
        x(i) = (1/A(i,i))*(b(i)-sigma);
       % x(i)= (1-w)*x(i) + (w)*(1/A(i,i))*(b(i)-sigma);
    end
    V = [V x];
    err = max(abs(x - x_old));
    residerr = max(abs(A*x - b))/norm(b);
    resid = [resid; residerr];
    relerr = [relerr; err];
    if residerr < tol
        break
    end   
end

sol = x;
P = V;
iterno = iter;
end



   


    