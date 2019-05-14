function [ parabolic ] = ParabolicCylinderFuns(ydata,modes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Create hermite polynomials in form of a matrix 
    
    Ank = zeros(modes);
    
    %first 2 polynomials are given, 
    Ank(1,1) = 1;
    Ank(1,2) = 0;
    Ank(2,2) = 2;
    
    %Now use recursion relation. Should give a lower triangular matrix. 
    for n = 3:modes
        %k = 1 case
        Ank(n,1) = -2*(n-2)*Ank(n-2,1);
        for k = 2:n
            Ank(n,k) = 2*Ank(n-1,k-1)-2*(n-2)*Ank(n-2,k);
        end
    end
    
    powersy = vander(ydata);
    powersy = fliplr(powersy);
    powersy = powersy(:,1:modes);
    
    parabolic = zeros(length(ydata),modes);
    
    for i = 1:modes
        coeffi = (factorial(i-1)*sqrt(pi))^(-1/2)*2^(-(i-1)/2);
        
        parabolic(:,i) = coeffi*sum((ones(length(ydata),1)*Ank(i,:)).*powersy,2).*exp(-ydata.^2/2);
    end
    


end

