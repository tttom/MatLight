% RGB=redHotColorMap(nbEntries)
%
% Creates a color map to be used with the command colormap
% All arguments are optional.
%
% Example usage:
%    figure;
%    image([0:.001:1]*256);
%    colormap(redHotColorMap(256))
function RGB=redHotColorMap(nbEntries)
    if (nargin<1 || isempty(nbEntries))
        nbEntries=64;
    end
    
    colors=[1 1 1; 1 1 0; 1 0 0;  0 0 0];
%     colorPositions=[0 .375 .75 1];
    colorPositions=[0 .375 .625 .95];
    RGB=interpolatedColorMap(nbEntries,colors,colorPositions);
end