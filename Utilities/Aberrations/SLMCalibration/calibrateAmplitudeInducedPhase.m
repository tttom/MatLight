%
% [graylevels,phases]=calibrateAmplitudeInducedPhase(slm,cam)
%
% This code is in a beta stage only!
%
% Attempts to measure the phase error when doing amplitude modulation on a Dual Head SLM.
%
% Input parameters:
%  slm: the DualHeadSLM object to measure
%  cam: the camera for the measurement
% 
% Output parameters:
%  graylevels: the measured graylevels, normalized to the dynamic range
%  phases: the measured phase error for each graylevel
%
function [graylevels,phases]=calibrateAmplitudeInducedPhase(slm,cam)
    if (nargin<1 || isempty(slm))
        slm=DualHeadSLM(2,[1/18 1/18]);
        slm.twoPiEquivalent=1.0;
    end
    referenceDeflectionFrequency=slm.referenceDeflectionFrequency;
    slm.referenceDeflectionFrequency=[0 0]; %Set to zero for now
    if (nargin<2 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=25e-3;
        cam.defaultNumberOfFramesToAverage=10;
    end
    
    roiSize=[128 128]; % [height, width]
            
    [X,Y]=meshgrid([1:slm.maxSize(2)],[1:slm.maxSize(1)]);
    referenceDeflection=exp(2i*pi*(X*slm.referenceDeflectionFrequency(2)+Y*slm.referenceDeflectionFrequency(1)));
                
    %Find the zeroth order spot center on the CCD
    slm.modulate(1);
    zerothOrderImage=cam.acquire();
    %Find the first order spot center on the CCD
    slm.referenceDeflectionFrequency=referenceDeflectionFrequency;
    slm.modulate(1);
    initialImage=cam.acquire()-zerothOrderImage;
    [ign, maxIdx]=max(initialImage(:));
    [maxRow, maxCol]=ind2sub(size(initialImage),maxIdx);
    
    centerPos=[maxRow,maxCol]; % [row, col]
    %Adjust the region of interest
    centerPos=min(max(centerPos,floor((roiSize+1)/2)),size(initialImage)-floor((roiSize+1)/2));%Make sure that the region of interest is on the CCD
    logMessage('Centering the region of interest around the coordinates [row,col]=[%u,%u].',centerPos);
        
    roi=[centerPos([2 1])-floor(roiSize([2 1])/2), roiSize([2 1])]; % [col, row, width, height]
    cam.regionOfInterest=roi;
    
    %Get dark image
    slm.modulate(0); %The CCD should be black in principle
    cam=cam.acquireBackground(100);
    
    graylevels=[0:16:255]./256;
    graylevels=fftshift(graylevels); %Start somewhere in the centre where there is signal
    values=zeros(size(graylevels));
    probeIdx=0;
    for (idx=1:length(graylevels))
        slm.modulate(graylevels(idx)*(angle(referenceDeflection)<0));
        img=cam.acquire();
%         values(idx)=max(img(:));
        if (probeIdx<=0)
            [ign probeIdx]=max(img(:));
        end
        valueTimesIntensity=median(img(probeIdx+[-size(img,1)-1 -size(img,1) -size(img,1)+1 -1 0 1 size(img,1)-1 size(img,1) size(img,1)+1]));
        
        %Compare with amplitude
        slm.modulate(graylevels(idx));
        img=cam.acquire();
        intensity=median(img(probeIdx+[-size(img,1)-1 -size(img,1) -size(img,1)+1 -1 0 1 size(img,1)-1 size(img,1) size(img,1)+1]));
        
        values(idx)=valueTimesIntensity./max(intensity,eps(1));
        
%         imagesc(img); colorbar()
%         drawnow;
    end
    graylevels=ifftshift(graylevels);
    values=ifftshift(values);
    values=max(0,values);%Negative values are just noise
    
    phases=acos(1-2*sqrt(values(beginIdx:endIdx)./max(values(:))));
    
    figure;
    plot(graylevels,phases);
    totalStroke=phases(end)+(phases(end)-phases(1))/(length(phases)-1) - phases(1);
    title(sprintf('total stroke: %.2fpi',totalStroke/pi));
    
end