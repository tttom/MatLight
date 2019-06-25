% [pupilFunctionWithoutDefocusTiltAndPiston defocus tilt piston]=removeDefocusTiltAndPiston(pupilFunction)
%
% Calculates and removes the defocus, tilt, and piston from a pupil function measurement.
%
% pupilFunctionWithoutTilt: the input with the tilt removed
% tilt: the removed tilt in units of Nyquist pixel shifts or number of waves over two whole aperture.
%       Divide this by size(pupilFunction) to get the shift in number of waves per pixel.
% piston: the removed piston in units of wavelength.
%
% Output parameters:
%     pupilFunctionWithoutTilt: the input with the tilt removed.
%     defocus: the removed defocus in number of wavelengths per pixel.
%     tilt: the removed tilt in units of wavelengths per pixel.
%     piston: the removed piston in units of wavelength.
%
% Example:
%     piston=0.25;
%     tilt=[0.05 -0.1];
%     defocus=0.004;
%
%     [X,Y]=ndgrid([-50:50],[-55:55]);
%     R2=X.^2+Y.^2;
%
%     pupilFunction=(sqrt(R2)<0.9*max(X(:))).*exp(2i*pi*piston).*exp(2i*pi*(tilt(1)*X-tilt(2)*Y + defocus*R2));
%     [pupilFunctionWithoutDefocusTiltAndPiston defocus tilt piston]=removeDefocusTiltAndPiston(pupilFunction);
%
%     figure();
%     subplot(1,2,1);
%     showImage(pupilFunction);
%     subplot(1,2,2);
%     showImage(pupilFunctionWithoutDefocusTiltAndPiston);
%
function [pupilFunctionWithoutDefocusTiltAndPiston defocus tilt piston]=removeDefocusTiltAndPiston(pupilFunction)
    if (nargin<1)
        % Test case
        [X,Y]=ndgrid([-50:50],[-55:55]);
        R2=X.^2+Y.^2;
        pupilFunction=exp(2i*pi*(0.05*X-0.1*Y + -.0004*R2))*1i.*(sqrt(R2)<0.5*max(X(:)));
    end
    
    inputSize=size(pupilFunction);
    
    weights=abs(pupilFunction);
    weightsX=weights(1:end-2,:).*weights(2:end-1,:).*weights(3:end,:);
    weightsY=weights(:,1:end-2).*weights(:,2:end-1).*weights(:,3:end);
    clear weights;
    RXX=[zeros(1,inputSize(2));exp(1i*diff(angle(pupilFunction),2,1)/4).*weightsX;zeros(1,inputSize(2))];
    RYY=[zeros(inputSize(1),1),exp(1i*diff(angle(pupilFunction),2,2)/4).*weightsY,zeros(inputSize(1),1)];
    RXXYY=RXX.*sign(RXX)+RYY.*sign(RYY); % angle mod pi
    meanAngle=angle(sum(RXXYY(:)));
    defocus=meanAngle/2/pi;
    
    [X,Y]=ndgrid([1:inputSize(1)]-1-floor(inputSize(1)/2),[1:inputSize(2)]-1-floor(inputSize(2)/2));
    R2=X.^2+Y.^2; clear X Y;
    pupilFunctionWithoutDefocus=pupilFunction.*exp(-2i*pi*R2*defocus);
    
    [pupilFunctionWithoutDefocusTiltAndPiston tilt piston]=removeTiltAndPiston(pupilFunctionWithoutDefocus);
end

%% Attempt to solve this in the Fourier domain
%     overSamplingRate=1;
%     
%     % linearize
%     maxR2=max(R2(:));
%     R2lin=[0:max(1,floor(sqrt(maxR2))):maxR2];
% %     R2lin=[0:maxR2];
%     Plin=zeros(1,length(R2lin)-1);
%     for RIdx=1:(length(R2lin)-1)
%         Plin(RIdx)=sum(pupilFunction(R2>=R2lin(RIdx) & R2<R2lin(RIdx+1)));
%     end
%     Plin(end)=Plin(end)+sum(pupilFunction(R2==R2lin(RIdx+1)));
%     R2lin=R2lin(1:end-1);
%     
%     Plin=Plin.*exp(-1i*angle(Plin(1))); % Remove phase offset at center
%     % oversample
%     if (overSamplingRate>1)
%         Plin(overSamplingRate*end)=0;
%         R2lin=diff(R2lin(1:2))*[0:(length(Plin)-1)];
%     end
%     
%     Plin=[Plin(1:end-1), conj(Plin(end:-1:2))];
%     R2lin=[R2lin(1:end-1), -R2lin(end:-1:2)];
%     
%     [ign I]=max(abs(ifft(Plin,'symmetric')));
%     defocus=-0.5/R2lin(I)/overSamplingRate;
% %     I=mod(I,length(Plin)/2);
%     
%     close all;
%     plot(fftshift(R2lin),fftshift(angle(Plin)), fftshift(R2lin),fftshift(abs(Plin))./max(abs(Plin(:))));
%     I    