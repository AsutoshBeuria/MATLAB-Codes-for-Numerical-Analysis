h =[0.1,0.2,0.3]; %step's size
for i=1
N =10;%number of steps
y = zeros(N+1, 1);
x = zeros(N+1, 1);
y(1)=1;
for i =1:N
    y(i+1) = y(i)+h*(-6*y(i));
    x(i+1) = i*h;
end
plot(x,y,'o')
hold on 
h=0.01; %step's size
N=10; %number of steps
y = zeros(N+1,1)
x = zeros(N+1,1) 
y(1) =1;
for i =1:N
   y(i+1) =  y(i)+h*(-6*y(i))
   x(i+1) = i*h;
end
plot(x,y)