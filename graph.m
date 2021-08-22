n=8;
x=zeros(n,1);
for i=1:n
    x(i)=10*i
end
for i=1:n
  y1(i)=2e-08*x(i)^5-3e-06*x(i)^4+0.0002*x(i)^3-0.0031*x(i)^2+0.0329*x(i)-0.0013
end

