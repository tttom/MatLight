% fieldValues=laguerreGaussian(R,P,Z,pValue,lValue,lambda,beamWaist)
%
% Input:
%     R: radial coordinate grid [m]
%     P: azimutal coordinate grid [rad]
%     Z: axial coordinate grid [m]
%     pValue: the p value of the beam (radial)
%     lValue: the l value of the beam (azimutal phase index, must be integer)
%     lambda: the wavelength to simulate [m]
%     beamWaist: the beam waist [m]
%
% Example:
%         xRange=[-4e-6:.05e-6:4e-6];
%         [X,Y]=meshgrid(xRange,xRange);
%         [P,R]=cart2pol(X,Y);
%         Z=[];
%         lambda=500e-9;
%         beamWaist=2*lambda;
%         pValue=3;
%         lValue=0;
%         fieldValues=laguerreGaussian(R,P,Z,pValue,lValue,lambda,beamWaist);
%         imagesc(abs(fieldValues).^2)
%
function fieldValues=laguerreGaussian(R,P,Z,pValue,lValue,lambda,beamWaist)
    if (nargin<1),
        xRange=[-4e-6:.05e-6:4e-6];
        [X,Y]=meshgrid(xRange,xRange);
        [P,R]=cart2pol(X,Y);
        Z=[];
        pValue=0;
        lValue=0;
    end

    if (nargin<6)
        lambda=500e-9;
    end
    if (nargin<7)
        beamWaist=2*lambda;
    end
    
    if (isempty(Z)),
        Z=zeros(size(R));
    end
    
    rayleighRange=pi*beamWaist^2/lambda;
    W=beamWaist*sqrt(1+(Z./rayleighRange).^2);
    radiusOfCurvature=(Z+(Z==0)).*(1+(rayleighRange./Z).^2);
    gouyPhase=atan2(Z,rayleighRange);
    fieldValues=(1./W).*(R.*sqrt(2)./W).^abs(lValue).*exp(-R.^2./W.^2)...
        .*L(abs(lValue),pValue,2*R.^2./W.^2).*exp(1i.*(2*pi/lambda).*R.^2./(2*radiusOfCurvature))...
        .*exp(1i*lValue.*P).*exp(-1i*(2*pValue+abs(lValue)+1).*gouyPhase);
    
    %Normalize
    fieldValues=fieldValues./sqrt(sum(abs(fieldValues(:)).^2));
    
    if (nargout==0)
        showImage(fieldValues+1i*sqrt(eps('single')),-1,X,Y); axis square;        
%         intensityValues=abs(fieldValues).^2;
%         beamWaistEstimate=sum(intensityValues(:).*R(:))./sum(intensityValues(:));        
%         psf=fftshift(abs(fft2(fieldValues)).^2); psf=exp(1)^2*psf/max(psf(:));
%         beamDisplacement=[1 0];
%         fieldValuesDisp=fieldValues.*exp(2i*(beamDisplacement(1)*X+beamDisplacement(2)*Y)/beamWaist);
%         psfDisp=fftshift(abs(fft2(fieldValuesDisp)).^2); psfDisp=exp(1)^2*psfDisp/max(psfDisp(:));
%         plot([psf(1+floor(end/2),:); psfDisp(1+floor(end/2),:)].')
        
        clear fieldValues;
    end
end

function Y=L(lValue,pValue,X)
    Y=polyval(LaguerreGen(pValue,lValue),X);
end