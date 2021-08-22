function d = eigenvals(A)
    maxit = 2;
    for i = 1:maxit
        [c,v,e,h,q,r] = QRfact(A)
        if i == 1
            c_01 = c;
            v_01 = v;
            e_01 = e;
            H_01 = h;
        end
        A = r*q
    end
    d = diag(A)
end