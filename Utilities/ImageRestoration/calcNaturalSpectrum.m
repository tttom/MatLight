% spectrum=calcNaturalSpectrum(xfRange,yfRange,kSignal)
%
% Calculate the spectrum amplitude at the frequencies on a regular plaid grid and averaging over each rectangular cell.
% The spectrum amplitude is defined as kSignal/frequency, where kSignal is
% a constant and frequency is magnitude of the spatial frequency.
% 
% Inputs:
%     x/yfRange: vectors with monotonously linearly increasing values indicating the
%     spatial frequencies of the 2D grid
%     kSignal: optional, default 1, signal strength multiplication factor
%
% Outputs:
%     spectrum: the average spectrum amplitude for each grid cell
%
function spectrum=calcNaturalSpectrum(xfRange,yfRange,kSignal)
    if nargin<3 || isempty(kSignal),
        kSignal=1.0;
    end
    
    cellSize=[diff(xfRange(1:2)) diff(yfRange(1:2))];
    
    xfRangeEdge=xfRange([1,1:end])+[-0.5, 0.5*ones(1,numel(xfRange))]*cellSize(1);
    yfRangeEdge=yfRange([1,1:end])+[-0.5, 0.5*ones(1,numel(yfRange))]*cellSize(2);
    [XFE,YFE]=ndgrid(xfRangeEdge,yfRangeEdge);
    [TFE,RFE]=cart2pol(XFE,YFE);
    spectrumIntegral=XFE.*log(RFE+YFE)+YFE.*log(RFE+XFE);

    spectrum=kSignal*diff(diff(spectrumIntegral,1,1),1,2)./prod(cellSize);
end