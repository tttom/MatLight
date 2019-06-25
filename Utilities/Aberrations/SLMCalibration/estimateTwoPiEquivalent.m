% [twoPiEquivalent slm]=estimateTwoPiEquivalent(slm,cam,probePos,probeSize)
%
% This code is in a beta stage only!
%
% Optimizes the mean intensity in the region of interest using direct search.
%
% Input parameters:
%  slm: the SLM object to measure
%  cam: the camera for the measurement
%  probePos: the position ([bottom,right]) of the focal spot where to measure the intensity change (default center)
%  probeSize: the area height and width to measure and take the median of (default [4 4]).
% 
% Output parameters:
%  twoPiEquivalent: the detected value indicating the fraction of the dynamic range that corresponds to 2pi modulation
%  slm: the input slm object with the twoPiEquivalent property updated
%
function [twoPiEquivalent slm]=estimateTwoPiEquivalent(slm,cam,probePos,probeSize)
    if (nargin<3 || isempty(probePos))
        probePos=cam.regionOfInterest(1:2)+1+floor(cam.regionOfInterest(3:4)/2);
    end
    if (nargin<4 || isempty(probeSize))
        probeSize=[1 1]*5;
    end
    
    minIntensity=0.01;
        
    testIntensity=1.0;
    while testIntensity>minIntensity && isAlmostSaturated(slm,cam,testIntensity)
        testIntensity=testIntensity/2;
    end
    
    if (testIntensity>minIntensity)
        [twoPiEquivalent,fval,exitflag,output] = fminsearch(@(x) -measureEfficiency(slm,cam,probePos,probeSize,testIntensity,x),slm.twoPiEquivalent,optimset('TolX',0.001,'MaxIter',20,'Display','iter'));

        slm.twoPiEquivalent=twoPiEquivalent;
    else
        twoPiEquivalent=slm.twoPiEquivalent;
        logMessage('Couldn''t get phase estimate due to image saturation!');
    end
end

function almostSaturated=isAlmostSaturated(slm,cam,testIntensity)
    slm.modulate(testIntensity);
    img=cam.acquire()+cam.background;
    almostSaturated=any(img(:)>.80);
end
    
function efficiency=measureEfficiency(slm,cam,probePos,probeSize,testIntensity,twoPiEquivalent)
    slm.twoPiEquivalent=twoPiEquivalent;
    slm.modulate(testIntensity);
    
    img=cam.acquire();
    probePos=probePos-cam.regionOfInterest(1:2);
    probeRange1=probePos(1)+[-floor(probeSize(1)/2):floor((probeSize(1)-1)/2)];
    probeRange2=probePos(2)+[-floor(probeSize(2)/2):floor((probeSize(2)-1)/2)];
    img=img(probeRange1(probeRange1>0 & probeRange1<=size(img,1)),probeRange2(probeRange2>0 & probeRange2<=size(img,2)),:);
    
    efficiency=median(img(:));
end

function value=defaultProbeFunctor(img)
    rng=[-5:5];
    
    if (size(img,1)>length(rng))
        img=img(1+floor(end/2)+rng,:,:);
    end
    if (size(img,2)>length(rng))
        img=img(:,1+floor(end/2)+rng,:);
    end
    
    value=median(img(:));
end