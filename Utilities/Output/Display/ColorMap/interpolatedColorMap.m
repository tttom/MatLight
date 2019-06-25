% RGB=interpolatedColorMap(nbEntries,colors,colorPositions)
%
% Creates a color map to be used with the command colormap
% All arguments are optional.
%
% Example usage:
%    figure;
%    image([0:.001:1]*256);
%    colormap(interpolatedColorMap(256,[0 .3 0; .7 .3 0; .7 .3 1; .7 1 1],[0 .3 .75 1]))
function RGB=interpolatedColorMap(nbEntries,colors,colorPositions)
    if (nargin<1 || isempty(nbEntries))
        nbEntries=64;
    end
    if (nargin<2 || isempty(colors))
        if (nargin<3 || isempty(colorPositions))
            colors=([1:length(colorPositions)]-1)/(length(colorPositions)-1);
        else
            colors=[0 0 0; 1 1 1];
        end
    end
    if (nargin<3 || isempty(colorPositions))
        colorPositions=[0:1/(size(colors,1)-1):1];
    end
    
    [colorPositions sortI]=sort(colorPositions);
    colors=colors(sortI,:);
    
    fraction=([1:nbEntries]-1)/max(1,nbEntries-1);
    
    RGB=interp1q(colorPositions.',colors,max(min(fraction.',colorPositions(end)),colorPositions(1)));
    
    RGB=min(1,max(0,RGB));
end