% coefficients=zernikeDecomposition(X,Y,Z,maxTerms)
%
% X, Y and Z should be matrices of real numbers.
%
% Returns:
%     coefficients: the vector of standard zernike coefficients, the first of
%                   which are: piston,
%                              tip(x),tilt(y),
%                              defocus, astigmatism-diag,astigmatism-x,
%                              coma-y,coma-x,  trefoil-y,trefoil-x,
%                              spherical aberration
%                   where the postscripts indicate the position of the extreme
%                   value on the pupil edge.
%
%
% See also: zernikeComposition.m, and zernike.m
%
function coefficients=zernikeDecomposition(X,Y,Z,maxTerms)
    R=sqrt(X.^2+Y.^2);
    T=atan2(Y,X);
    invalidPoints=any(isnan(Z(R<=1)));
    for j=1:maxTerms,
        currZernike=zernike(j,R,T);
        if (invalidPoints || j==1),
            normalization=sum(currZernike(~isnan(Z)).^2);
        end
        coefficients(j)=sum(Z(~isnan(Z)).*currZernike(~isnan(Z)))./normalization;
    end
end