function [L,U,P] = lu_dcmp(A)
%This gives LU decomposition of A with the permutation matrix P
% denoting the row switch(exchange) during factorization
NA = size(A,1);
AP = [A eye(NA)]; %augment with the permutation matrix.
for k = 1:NA - 1
%Partial Pivoting at AP(k,k)
[akx, kx] = max(abs(AP(k:NA,k)));
if akx < eps
error('Singular matrix and No LU decomposition')
end
mx = k+kx-1;
if kx > 1 % Row change if necessary
tmp_row = AP(k,:);
AP(k,:) = AP(mx,:);
AP(mx,:) = tmp_row;
end
% LU decomposition
for m = k + 1: NA
AP(m,k) = AP(m,k)/AP(k,k); %Eq.(2.4.8.2)
AP(m,k+1:NA) = AP(m,k + 1:NA)-AP(m,k)*AP(k,k + 1:NA); %Eq.(2.4.9)
end
end
P = AP(1:NA, NA + 1:NA + NA); %Permutation matrix
for m = 1:NA
for n = 1:NA
    if m == n, L(m,m) = 1.; U(m,m) = AP(m,m);
elseif m > n, L(m,n) = AP(m,n); U(m,n) = 0.;
else L(m,n) = 0.; U(m,n) = AP(m,n);
end
end
end
if nargout == 0, disp('L*U = P*A with'); L,U,P, end
%You can check if P?*L*U = A?
