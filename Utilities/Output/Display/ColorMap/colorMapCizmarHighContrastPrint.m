% RGB=colorMapCizmarHighContrastPrint(nbEntries)
%
% Based on Tomas' colormap for producing good contrast black and white printouts.
% Creates a color map to be used with the command colormap.
% All arguments are optional.
%
% Example usage:
%    figure;
%    image([0:.001:1]*256);
%    colormap(colorMapCizmarHighContrastPrint(256))
function RGB = colorMapCizmarHighContrastPrint(nbEntries)
    if (nargin<1 || isempty(nbEntries))
        nbEntries=64;
    end
    
    RGB=interpolatedColorMap(nbEntries,[0 0 .3; .7 0 .3; .7 1 .3; .7 1 1],[0 .292 .708 1]);
end

% levs=levs(:,1);
% t=[min(2.4*levs,.7),levs,max(2.4*levs-1.4,.3)];
% t(:,2)=2.4*t(:,2)+.3-t(:,1)-t(:,3);
% t=max(0,min(1,t));