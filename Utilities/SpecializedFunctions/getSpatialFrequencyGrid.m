% varargout=getSpatialFrequencyGrid(pixelPitch,gridSize)
%
% Creates N matrices of size gridSize with the spatial frequency coordinates for a
% sample grid with pixel separation pixelPitch. gridSize and pixelPitch must be vectors of equal length N.
%
% Examples:
%    [XOtf YOtf]=getSpatialFrequencyGrid([1e-6 1e-6],[512 512]);    
%    cutOffSpatialFrequency=(2*config.detection.objective.numericalAperture)/config.detection.wavelength;
%    fRel=sqrt(XOtf.^2+YOtf.^2)./cutOffSpatialFrequency;
%
%    [XOtf YOtf ZOtf]=getSpatialFrequencyGrid([1 1 1]*100e-9,[1 1 1]*128);
%
function varargout=getSpatialFrequencyGrid(pixelPitch,gridSize)
    nbDims=numel(pixelPitch);
    sampleFrequencies=1./pixelPitch;
    otfSteps=sampleFrequencies./gridSize;
    freqRange=cell(1,nbDims);
    for dimIdx=1:nbDims,
        freqRange{dimIdx}=([0:gridSize(dimIdx)-1]-floor(gridSize(dimIdx)/2))*otfSteps(dimIdx);
    end
    
    freqRange{end+1}=[0];
    
    varargout=cell(1,nargout);
    [varargout{:}]=ndgrid(freqRange{:});
end