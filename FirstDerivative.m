function [deriv1] = FirstDerivative(X,Y)
% INPUTS:
% X and Y are vectors containing the coordinates of given data points 
%OUTPUT:
% deriv1 :vector of first derivative of the function described by the data set, at ALL given data points

deriv1=zeros(length(X),1);
n = length(X);
for i = 1:n-2
    for j = 1:n-2
        t = (2*X(i) - X(i+1) - X(i+2))*Y(i)/((X(i) - X(i+1)*(X(i) - X(i+2)) + (X(i) - X(i+2))*Y(i+1)/((X(i+1) - X(i)*(X(i+1) - X(i+2)) + (X(i) - X(i+1))*Y(i+2)/((X(i+2) - X(i)*(X(i+2) - X(i+1));
        deriv1 = [deriv1 ; t];
    end
    
end




end

