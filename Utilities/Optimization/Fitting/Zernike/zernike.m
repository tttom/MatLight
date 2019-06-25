%
% Returns samples of a sqrt(pi)*orthonormal basis on the unit disk,
% for 2x2-unit square, multiply with 4/pi.
%
% WARNING: returns two complementary zernike polynomials at the same time as a complex function!
% User real(zernike(...)) to only get the requested one.
%
% result=zernike(n,m,rho,theta)
%    For m>=0, returns the even Zernike polynomial(cos) value as the real part,
%    and the odd polynomial(sin) value as the imaginary part. For m<0, the odd
%    Zernike value is returned as the real part, and the even is returned
%    as the imaginary part.
%    prerequisite: rho>=0
%
% result=zernike(j,rho,theta)
%    Returns the Zernike polynomial with standard coefficient j
%       "Zernike polynomials and atmospheric turbulence", Robert J. Noll,
%       JOSA, Vol. 66, Issue 3, pp. 207-211, doi:10.1364/JOSA.66.000207
%     The first of which are: piston,
%                             tip(x),tilt(y),
%                             defocus, astigmatism-diag,astigmatism-x,
%                             coma-y,coma-x,  trefoil-y,trefoil-x,
%                             spherical aberration
%                  where the postscripts indicate the position of the extreme
%                  value on the pupil edge.
%
%This function can handle rho and theta matrices and n and m vectors.
%
%
% Example:
%     gridSize=[1 1]*512;
%     uRange=-1:(2/(gridSize(2)-1)):1;
%     vRange=-1:(2/(gridSize(1)-1)):1;
%     [X,Y]=meshgrid(uRange,vRange);
%     R=sqrt(X.^2+Y.^2);
%     T=atan2(Y,X);
%     figure;
%     ssurf(real(zernike(4,R,T))./(R<1));
%     figure;
%     subplot(1,2,1);
%     ssurf(real(zernike(3,1,R,T))./(R<1));
%     subplot(1,2,2);
%     ssurf(imag(zernike(3,1,R,T))./(R<1));
%
%     mean(mean(real(zernike(3,1,R,T)).*real(zernike(2,0,R,T)).*(R<1)))*4/pi
%     mean(mean(real(zernike(3,1,R,T)).*real(zernike(3,1,R,T)).*(R<1)))*4/pi
%
% See also: zernikeFit.m and zernikeComposition.m as well as utility
% functions ind2subZernike.m and sub2indZernike.m
%
function result=zernike(n,m,rho,theta)
    wantStandardCoefficientOnly=nargin<4;
    if (wantStandardCoefficientOnly),
        if nargin<3,
            theta=[];
        else
            theta=rho;
        end
        rho=m;
        js=n;%Standard Zernike coefficient number:
        [n,m]=ind2subZernike(js);
    end
    
    dataSize=size(rho);
    nbModes=numel(m);
    
    if isempty(theta),
        theta=zeros(dataSize);
    end
    
    normalization=sqrt(2*(n+1)./(1+(m==0)));%Make orthonormal basis on unit disk (for 2x2-unit square, multiply with 4/pi)
    result=zeros([dataSize nbModes]);
    for modeIdx=1:numel(m),
        if mod(n-m,2)==0,
            zernikePhi=exp(1i*(m(modeIdx)*theta+(m(modeIdx)<0)*pi/2)); % Set the real part as requested
            result([1:numel(rho)]+(modeIdx-1)*numel(rho))=normalization(modeIdx)*zernikeR(abs(m(modeIdx)),n(modeIdx),rho).*zernikePhi;
        end
    end
end

% Determine the radial component, the zernike polynomial, for all rho in a matrix
% prerequisites: m>=0, rho>=0, mod(n-m,2)==0
% Output: a matrix of size rho, or the scalar 0 indicating an all zero
% result in case the difference n-m is odd.
%
function result=zernikeR(m,n,rho)
    result=0;
%     rho=(rho<=1).*rho;
    if (mod(n-m,2)==0),
        rhoPow=rho.^m;
        rhoSqd=rho.^2;
        
        lCoeffs=[((n-m)/2):-1:0];
        % largest factorial result is n!
        % subResultWeights=(((-1).^lCoeffs).*factorial(n-lCoeffs))./(factorial(lCoeffs).*factorial((n+m)/2-lCoeffs).*factorial((n-m)/2-lCoeffs));
        
        for lIdx=1:numel(lCoeffs),
            lCoeff=lCoeffs(lIdx);
            %For speedup: rhoPow=rho.^(n-2*lCoeffs);
            if lIdx>1,
                rhoPow=rhoPow.*rhoSqd; % note the lCoeffs are in reversed order
            end
            subResultWeight=((-1)^lCoeff).*factorialFraction(n-lCoeff,[lCoeff (n+m)/2-lCoeff (n-m)/2-lCoeff]);
            result=result+subResultWeight*rhoPow;
        end
    end
end