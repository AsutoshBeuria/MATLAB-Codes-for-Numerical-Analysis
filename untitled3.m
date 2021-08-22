h =0.1; %step's size
N =10;%number of steps
y = zeros(N+1, 1);
x = zeros(N+1, 1);
y(1)=1;
for i =1:N
    y(i+1) = y(i)+h*(-6*y(i));
    x(i+1) = i*h;
end
plot(x,y)
h=0.01; %step's size
N=10; %number of steps
y = zeros(N+1,1)
x = zeros(N+1,1) 
y(1) =1;
for i =1:N
   y(i+1) = y(i)+h*(-6*y(n))
   x(i+1) = n*h;
end
plot(x,y)