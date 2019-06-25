%The real part of the coefficients are for the even polynomials,
%the imaginary part of the coefficients are for the odd polynomials
function result=zernikeLinearCombination(coefficients,X,Y)
    result=zeros(size(X));
    [M,N]=getMNFromCoeffiecientUpToIndex(length(coefficients));
    for index=1:length(coefficients)
        if (abs(coefficients(index))>0)
            subResult=zernike(M(index),N(index),sqrt(X.^2+Y.^2),atan2(Y,X));
            result=result + real(subResult)*real(coefficients(index)) + imag(subResult)*imag(coefficients(index));
        end
    end
end

% n>=m>=0 and n-m even
% Not the official order, but:
% (m,n) = (0,0) (1,1) (0,2) (1,1) (1,3) (2,2) (0,4) (1,3) (2,2) (1,5) (2,4)
%          (3,3) (0,6) (1,5) (2,4) (3,3) (1,7) (2,6) (3,5) (4,4)
function [M,N]=getMNFromCoeffiecientUpToIndex(coeffIndex)
    M=[];
    N=[];
    l=0;
    m=0;
    n=0;
    for index=1:coeffIndex-1,
        M(end+1)=m;
        N(end+1)=n;
        if (m<n)
            m=m+1;
            n=n-1;
        else
            l=l+1;
            n=l;
            m=mod(n,2);
        end
    end
    M(end+1)=m;
    N(end+1)=n;
end