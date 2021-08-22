function x = Gauss(A,b)
% create the augmented matrix
Ab = [A,b];
% dtermine the size of the augmented matrix
[R,C] = size(Ab);

% iterate columns
for i=1:R-1
    
    % pivoting
    if(abs(Ab(i,i)) ~= max(abs(Ab(i:R,i))))
        for j = i+1:R
            if(abs(Ab(j,i)) == max(abs(Ab(i:R,i))))
                temp = Ab(i,:) ;
                Ab(i,:) = Ab(j,:) ;
                Ab(j,:) = temp ;
                break;
            end
            x = Ab
        end
    end
 
    % Gaussian elimination
    for k = i+1:R
        Ab(k,i:C) = Ab(k,i:C) - Ab(k,i)/Ab(i,i)*Ab(i,i:C);
    end 
    y = Ab
end

% take the the upper triangular matrix and vector out from Ab
U = Ab(1:R,1:R);
b = Ab(:, R+1);

% Use the back substitution to solve the remaining equations
x = backSub(U, b);
% if you have trouble with creating the function backSub(), you can
% romove line 36 and activate the line below.
% x = U \ b;
end