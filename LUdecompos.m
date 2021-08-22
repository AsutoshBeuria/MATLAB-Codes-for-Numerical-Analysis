function [L,U] = LUdecompos(A,b)
% create the augmented matrix
Ab = [A,b];
% dtermine the size of the augmented matrix
[R,C] = size(Ab);
L = eye(R);
% iterate columns
for i=1:R-1
    
    % pivoting
    if(abs(Ab(i,i)) ~= max(abs(Ab(i:R,i))))
        for j = i+1:R
            if(abs(Ab(j,i)) == max(abs(Ab(i:R,i))))
                temp = Ab(i,:) ;
                Ab(i,:) = Ab(j,:) ;
                Ab(j,:) = temp ;
                Ab
                L([i,j],1:i-1) = L([j,i],1:i-1);
                L
                break;
            end
        end
    end
 
    % Gaussian elimination
    for k = i+1:R
        L(k,i)= Ab(k,i)/Ab(i,i)
        Ab(k,i:C) = Ab(k,i:C) - L(k,i)*Ab(i,i:C)
    end 
end

% take the the upper triangular matrix and vector out from Ab
U = Ab(1:R,1:R)

end