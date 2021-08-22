function no = Jacobi(A,b,imax,x0,tol)
x = zeros(imax,length(b));

% Iteration k = 1
x(1,:) = x0';

for kk = 2:imax
    for ii = 1:length(b)
        vec = 1:length(b);
        vec(vec==ii) = [];
        x(kk,ii) = (1/A(ii,ii))*(b(ii) - sum(A(ii,vec).*x(kk-1,vec)));
    end
    error_metric(kk-1,:) = norm( x(kk,:) -x(kk-1,:) );
    if error_metric(kk-1,:) < tol
        break
    end
end
X = x(1:kk,:);
array2table(X);
x_Jacobi = X(kk,:).';
no = kk;
end
