% [cam centerPos]=selectRegionOfInterestAroundPeakIntensity(cam,slm,roiSizeForCam,centerPos)
%
% Updates the cameras cam's region of interest to a size roiSizeForCam (rows, columns) around the peak intensity, and re-acquires the dark image. 
%
% Argument centerPos is optional, when given the algorithm assumes that the
% peak intensity is at the coordinates as given by centerPos.
function [cam centerPos]=selectRegionOfInterestAroundPeakIntensity(cam,slm,roiSizeForCam,centerPos)
    if (nargin<4 || isempty(centerPos))
        %Capture dark image
        slm.modulate(0); %The CCD should now be unilluminated in principle
        cam=cam.acquireBackground(cam.numberOfFramesToAverage*4);
        
        irradiationIntensityScaling=1;
        maxValue=2;
        while(maxValue>.90 && irradiationIntensityScaling>1/16)
            %Capture the spot without the zeroth order
            slm.modulate(irradiationIntensityScaling);
            initialImage=cam.acquire();

            %Get peak position
            [maxValue, maxIdx]=max(initialImage(:));
            
            irradiationIntensityScaling=irradiationIntensityScaling/2;
        end
        [maxRow, maxCol]=ind2sub(size(initialImage),maxIdx);
        centerPos=[maxRow,maxCol]; % [row, col]
        centerPos=centerPos+cam.regionOfInterest(1:2); %There might be an initial offset
        
        logMessage('Centering the camera region of interest around the coordinates [row,col]=[%u,%u].',centerPos);
    else
        logMessage('Forcing the camera region of interest to be around the coordinates [row,col]=[%u,%u].',centerPos);
    end
    
    %Adjust the region of interest
    centerPos=2*floor(centerPos/2); %Align with Bayer filter if necessary
    %Make sure that the region of interest is on the CCD
    centerPos=min(max(centerPos,floor((roiSizeForCam+1)/2)),cam.maxSize-floor((roiSizeForCam+1)/2));
    cam.regionOfInterest=[centerPos-floor(roiSizeForCam/2), roiSizeForCam];
    slm.modulate(0);
    cam=cam.acquireBackground(cam.numberOfFramesToAverage*4);%Measure background for the new ROI
end