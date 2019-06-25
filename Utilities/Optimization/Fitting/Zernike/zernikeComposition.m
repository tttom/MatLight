% Z = zernikeComposition(X,Y,coefficients)
%
% X and Y should be matrices of real numbers.
% coefficients: the row vector of standard zernike coefficients, the first of
%               which are: piston,
%                          tip(x),tilt(y),
%                          defocus, astigmatism-diag,astigmatism-x,
%                          coma-y,coma-x,  trefoil-y,trefoil-x,
%                          spherical aberration
%               where the postscripts indicate the position of the extreme
%               value on the pupil edge.
%               When a matrix of multiple rows is supplied, a 3D stack of
%               Zernike compositions is return, one 2D matrix for each row.
%
% See also: zernikeFit.m, and zernike.m
%
function Z = zernikeComposition(X,Y,coefficients)
    nbCoefficients = size(coefficients,1);
    coefficients = permute(coefficients, [3 2 1]);
    Z=zeros([size(X), nbCoefficients]);
    [T,R]=cart2pol(X,Y);
    for j=1:size(coefficients,2),
        if any(coefficients(1,j,:) ~= 0),
            Z = Z + bsxfun(@times, coefficients(1,j,:), real(zernike(j,R,T)));
        end
    end
end