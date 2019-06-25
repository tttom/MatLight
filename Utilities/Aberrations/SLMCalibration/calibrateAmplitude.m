%
% [graylevels,amplitudes]=calibrateAmplitude(slm,cam)
%
% This code is in a beta stage only!
%
% Attempts to measure the amplitude response for a Dual Head SLM.
%
% Input parameters:
%  slm: the DualHeadSLM object to measure
%  cam: the camera for the measurement
% 
% Output parameters:
%  graylevels: the measured graylevels, normalized to the dynamic range
%  amplitudes: the measured amplitude for each graylevel (square root of the intensity)
%
function [graylevels,amplitudes]=calibrateAmplitude(slm,cam)
    if (nargin<1 || isempty(slm))
        slm=DualHeadSLM(2,[1/18 1/18]);
    end
    referenceDeflectionFrequency=slm.referenceDeflectionFrequency;
    slm.referenceDeflectionFrequency=[0 0]; %Set to zero for now
    if (nargin<2 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=25e-3;
        cam.defaultNumberOfFramesToAverage=10;
    end

    origTwoPiEquivalent=slm.twoPiEquivalent;
    slm.twoPiEquivalent=1.0;
    
    roiSize=[128 128]; % [height, width]
                 
    %Find the zeroth order spot center on the CCD
    slm.modulate(1);
    zerothOrderImage=cam.acquire();
    %Find the first order spot center on the CCD
    slm.referenceDeflectionFrequency=referenceDeflectionFrequency; %Set the deflection frequency we are interested in
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
    
    graylevels=[0:4:255]./256;
    graylevels=fftshift(graylevels); %Start somewhere in the center where there is signal
    intensities=zeros(size(graylevels));
    probeIdx=0;
    for (idx=1:length(graylevels))
        slm.modulate(graylevels(idx));
        img=cam.acquire();
%         intensities(idx)=max(img(:));
        if (probeIdx<=0)
            [ign probeIdx]=max(img(:));
        end
        intensities(idx)=median(img(probeIdx+[-size(img,1)-1 -size(img,1) -size(img,1)+1 -1 0 1 size(img,1)-1 size(img,1) size(img,1)+1]));
        
%         imagesc(img); colorbar()
%         drawnow;
    end
    slm.twoPiEquivalent=origTwoPiEquivalent;

    graylevels=ifftshift(graylevels);
    intensities=ifftshift(intensities);
    intensities=max(0,intensities);%Negative intensities are just noise
    
    amplitudes=sqrt(intensities);
    
    plot(graylevels,amplitudes);
    
end

