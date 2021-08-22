 function [c_0,v_0,e_0,H_0,Q,R] = QRfact(A)
    n = size(A,1);
    a = A;
    c = a(:,1);
    e = zeros(n,1);
    I = eye(n); 
    if c(1) > 0
        e(1) = 1;
    else
        e(1) = -1;
    end
    v_0 = c + norm(c,2)*e;
    c_0 = c; e_0 = e;
    H_0 = I - 2*(v_0*v_0')/(v_0'*v_0);
    Q = H_0; R = H_0*a;
    for i = 2:n-1
        e = zeros(n,1);
        c = R(:,i);
        c(1:i-1) = 0;
        if c(i) > 0
            e(i) = 1;
        else
            e(i) = -1;
        end
        v = c + norm(c,2)*e
        H = I - 2*(v*v')/(v'*v)
        Q = Q*H
        R = H*R
    end
    
end