function [sol,P_GS,iterno] = GaussModified(A,b,n_max,x0,tol) 
% Iteration k = 1 ;
x = x0; 
n = size(A,1);
P_GS = eye(n);
Ab =[A,b];
for i = 1:n-1
   if(Ab(i,i) == 0)
        for j = i+1:n
            temp = P_GS(i,:);
            P_GS(i,:) = P_GS(j,:);
            P_GS(j,:) = temp;
            temp = Ab(i,:) ;
            Ab(i,:) = Ab(j,:) ;
            Ab(j,:) = temp ;
            break;
        end
   end   
end

B = Ab(:,1:n);
b = Ab(:,end);
for iter = 1:n_max 
    x_old = x;
    for i = 1:n
        sigma = 0;      
        for j = 1:i-1
            sigma=sigma+B(i,j)*x(j);
        end
        for j=i+1:n
                sigma=sigma+B(i,j)*x_old(j);
        end     
        x(i) = (1/B(i,i))*(b(i)-sigma);
        
    end
    
    residerr = max(abs(B*x - b))/norm(b);
    if residerr < tol
        break
    end   
end

sol = x;
iterno = iter;
end
