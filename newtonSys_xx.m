function [v,k,x] = newtonSys_xx(J, fun, x0, xtol, ftol, maxit)
% newtonSys  Newton's method for systems of nonlinear equations.
%
% Synopsis:  x = newtonSys(J,fun,x0)
%            x = newtonSys(J,fun,x0,xtol)
%            x = newtonSys(J,fun,x0,xtol,ftol)
%            x = newtonSys(J,fun,x0,xtol,ftol,verbose)
%            x = newtonSys(J,fun,x0,xtol,ftol,verbose,arg1,arg2,...)
%
% Input:  J    = (string) name of mfile that returns the Jacobian matrix J
%         fun  = (string) name of mfile that returns vector f
%         x0   = initial guess at solution vector, x
%         xtol = (optional) tolerance on norm(dx).
%         ftol = (optional) tolerance on norm(f).
%
% Output:  x = solution vector;  x is returned after k iterations if tolerances
%              are met, or after maxit iterations if tolerances are not met.

% no
f_vec = feval(fun,x0);
if length(x0) ~= length(f_vec)
    disp('Error: length of the vector of initial estimates has to be: ')
    length(f_vec)
    return
end
x = x0;  % Initial guess and current number of iterations
v = [x0'];
for i = 1 : maxit
    k = i;
    J_mat = feval(J,x);   % Returns Jacobian matrix
    f_vec = feval(fun,x); % Returns f vector
    del_x = linsolve(J_mat,-f_vec);
    xipl = x + del_x;
    %x_err = max(abs(xipl - x));
    x_err = norm(del_x);
    f_err = abs(norm(feval(fun,xipl)));
   
    if x_err < xtol || f_err < ftol
        break
    end
   
    x = xipl;
    v = [v ; x'];
    
end
    %% The crucial part maybe lost after this line
    % However, to record the data or declare other variables you want, you 
    % are freely to add codes before this line.
    % --------------------------------------------------------------------------
    

end