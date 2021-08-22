function x = backSub(U, b)
n = size(U,1);

x(n,1) = b(n,1)/U(n,n);
for i = n-1:-1:1
    h = 0;
    for j = i:n-1
        h = h + U(i,j+1)*x(j+1,1);
    end
    
    x(i,1) = (b(i,1) - h)/U(i,i) ;
end

end