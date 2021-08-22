function x = newtonSys_modified(A, fun, x0, xtol, ftol, maxit)

x = x0;  % Initial guess and current number of iterations

for i = 1 : maxit
    J_mat = A;   % Returns Jacobian matrix
    n = cond(J_mat,2);
    if n > 1e4
        disp('Warning , the condition number is above 1e4')
    end
    
    if n > 1e10
        disp('The conditon no is above 1e10')
        break
    end
    
    f_vec = fun; % Returns f vector
    del_x = linsolve(J_mat,-f_vec);
    xipl = x + del_x;
    %x_err = max(abs(xipl - x));
    x_err = norm(del_x);
    f_err = abs(norm(fun));
   
    if x_err < xtol || f_err < ftol
        break
    end
   
    x = xipl;
end

end