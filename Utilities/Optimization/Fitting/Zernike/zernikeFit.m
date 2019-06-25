% [coefficients, residue, err]=zernikeFit(X,Y,Z,nbCoefficients,radialOnly,weights)
%
% Fits the first nbCoefficients Zernike polynomials to Z, defined in X, Y.
%
% X, Y and Z should be matrices of real numbers.
%
% Returns:
%     coefficients: the vector of standard Zernike coefficients, the first of
%                   which are: piston,
%                              tip(x),tilt(y),
%                              defocus, astigmatism-diag,astigmatism-x,
%                              coma-y,coma-x,  trefoil-y,trefoil-x,
%                              spherical aberration
%                   where the postscripts indicate the position of the extreme
%                   value on the pupil edge.
%
%
% See also: zernikeComposition.m, zernikeDecomposition.m, and zernike.m
%
function [coefficients, residue, err]=zernikeFit(X,Y,Z,nbCoefficients,radialOnly,weights)
    if nargin<2,
        gridSize=[1 1]*128;
        xRange=-1:(2/(gridSize(2)-1)):1;
        yRange=-1:(2/(gridSize(1)-1)):1;
        [X,Y]=meshgrid(xRange,yRange);
    end
    if nargin<3,
        Z=sqrt(X.^2+Y.^2);
    end
    if nargin<4 || isempty(nbCoefficients),
        nbCoefficients=100;
    end
    if nargin<5 || isempty(radialOnly),
        radialOnly=false;
    end
    if nargin<6 || isempty(weights),
        weights=ones(size(Z));
    end
    
    [T,R]=cart2pol(X,Y);
    inUnitDisk = R<1;
    validPos = ~isnan(Z) & inUnitDisk & weights>0;
    
    if ~radialOnly,
        zernikeI=[1:nbCoefficients];
    else
        nbRadialTerms=ind2subZernike(nbCoefficients)/2+1;
        zernikeI=sub2indZernike(2*([1:nbRadialTerms]-1));
    end
    zernikeBasis=real(zernike(zernikeI,R(validPos),T(validPos)));
    zernikeBasis=zernikeBasis(:,:);
    normalization=sqrt(numel(validPos).*pi/4); %sqrt(sum(inUnitDisk(:))); % orthonormalize (bar the invalid positions)
    zernikeBasis=zernikeBasis./normalization; % so that these square-sum to 1
    
    % Calculate the matrix for the weighted problem
    zernikeBasisWeighted=zernikeBasis.*repmat(weights(validPos),[1 numel(zernikeI)]);
    
    coefficients=zeros(1,nbCoefficients);
    [coefficients(zernikeI) flag]=minres(zernikeBasisWeighted'*zernikeBasisWeighted,zernikeBasisWeighted'*(weights(validPos).*Z(validPos)./normalization),1e-6,100);
    %err=zernikeBasisWeighted*coefficients(zernikeI).'-(Z(validPos)./normalization);  norm(err)
    % Minimize |A*Ax-A*B|^2 + r2|x|^2 => H=A*A: (H^2+r2)x == HA*B
%     H=zernikeBasis'*zernikeBasis;
%     r2=0;
%     randn('seed',0);
%     Zm=Z+.1*randn(size(Z));
%     coefficients(zernikeI)=minres(@(x) zernikeBasis'*zernikeBasis*x,zernikeBasis'*Z(validPos));
    %coefficients(zernikeI)=minres(H'*H+r2,H'*zernikeBasis'*Z(validPos));
    
    if nargout>=2,
        residue=Z(:)-squeeze(real(zernike(zernikeI,R(:),T(:))))*coefficients(zernikeI).';
        residue=reshape(residue,size(Z));
    end
    if nargout>=3,
        err=sqrt(sum(abs((residue(validPos)./normalization).*weights(validPos)).^2)./sum(abs(weights(validPos)).^2));
    end
    
    if nargout<1,
        zernikeFitted=zernikeComposition(X,Y,coefficients);
        
        figure('Name',sprintf('%0.3f ',coefficients));
        axs(1)=subplot(1,3,1);
        ssurf(X,Y,Z); title('surface to fit')
        axs(2)=subplot(1,3,2);
        ssurf(X,Y,zernikeFitted.*validPos./validPos); title('fitted');
        axs(3)=subplot(1,3,3);
        ssurf(X,Y,zernikeFitted.*validPos./validPos-Z);
        msg=sprintf('fit error on unit disk: %0.3f%%',100*sqrt(mean(abs(zernikeFitted(validPos)-Z(validPos)).^2)));
        title(msg);
        logMessage(msg);
        
        clear coefficients;
    end
end


