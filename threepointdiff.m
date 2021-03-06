function diff = threepointdiff(delx,f)
    y = f;
    n = size(y,1);
    d = zeros(n,1);
    for i = 1:n-2
        d(i) = (y(i+2) - 2*y(i+1) + y(i))/(delx^2);
    end
    diff = d;
end