clear all
close all
maxit = 1e3;
n = 6;
M = 0.5*eye(n);
K = zeros(n);
K(1:1+n:n*n) = 12e4;
K(n+1:1+n:n*n) = -6e4;
K(2:1+n:n*n-n) = -6e4;
K(n,n) = 6e4;
A = M\K;
for i = 1:maxit
    [c,v,e,h,q,r] = QRfact(A);
    if i == 1
        c_01 = c
        v_01 = v
        e_01 = e
        H_01 = h
    end
    A = r*q;
    if max(abs(tril(A,-1))) < 0.1
        c_0k = c
        v_0k = v
        e_0k = e
        H_0k = h
        Q = q
        R = r
        break
    end
end
no_of_iter = i
eigen_vals = diag(A)

%%%%%%%%%%%%%%%%%% Q2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K2 = K;
K2(n-1,n-1) = 102000; K2(n-1,n) = -42e3; K2(n,n-1) = -42e3; K2(n,n) = 42e3;
A = M\K;
A2 = M\K2;
[Phi, lambda] = eig(A);
[Phi2, lambda2] = eig(A2);

lambda = diag(lambda)
lambda2 = diag(lambda2)
phi = Phi
phi2 = Phi2

f = sqrt(lambda)./(2*pi)
f2 = sqrt(lambda2)./(2*pi)

x = linspace(1,n,n);
for i = 1:n
    figure();
    plot(x,phi(:,i),'-o',x,phi2(:,i),'-o');
    grid
    title("Mode shape",i)
    xlabel('Degree of freedom'); ylabel('Amplitude');
    legend('Phi(i)','Phi2(i)', 'Location','NorthOutside');
end









