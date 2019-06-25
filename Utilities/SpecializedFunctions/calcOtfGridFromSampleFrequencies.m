% [XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(sampleFrequencies,gridSize,cutOffSpatialFrequency)
%
% sampleFrequencies can be specified as a scalar for both dimensions or as a vector in the order [x y].
% gridSize must thus be specified in the order [y x] instead!
% 
% The size of the output arguments equals gridSize.
%
% Example:
%    pixelPitch=[1e-6 1e-6];
%    [XOtf,YOtf]=calcOtfGridFromSampleFrequencies(1./pixelPitch,[200 200]);
% or:
%    cutOffSpatialFrequency=5e5;
%    [XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(1./pixelPitch,[200 200],cutOffSpatialFrequency);
function [XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(sampleFrequencies,gridSize,cutOffSpatialFrequency)
    logMessage('Using obsolete function calcOtfGridFromSampleFrequencies.m. This will break in future releases. Please inform Tom Vettenburg tv2@st-andrews.ac.uk ');
    if (length(gridSize)==1),
        gridSize(2)=gridSize(1);
    end
    otfSteps=sampleFrequencies./gridSize([2 1]); %In order [x y]
    rangeX=([0:gridSize(2)-1]-floor(gridSize(2)/2))*otfSteps(1);
    rangeY=([0:gridSize(1)-1]-floor(gridSize(1)/2))*otfSteps(2);
    [XOtf,YOtf]=meshgrid(rangeX,rangeY);
    if (nargin>2),
        fRel=sqrt(XOtf.^2+YOtf.^2)./cutOffSpatialFrequency;
    end
end