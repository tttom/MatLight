% values=falseColorImageToValues(imgRGB,colorMap,range,defaultValue)
%
% imgRGB: a color image with the same number as color channels as
% represented by colorMap
% colorMap: a matrix with the color index (row idx) and channel value (col
% idx).
% range: optional vector of start and end value, or one value per colormap entry.
% defaultValue: optional default value if no nearby color is found, default NaN
%
% Converts a false color image and a colormap back to the underlying
% values. It uses the 2-norm on the RGB colors to find the closest match.
% Colors that are further separated than the largest difference in the
% colormap are set to the default value.
%
function values=falseColorImageToValues(imgRGB,colorMap,range,defaultValue)
    imgSize=size(imgRGB);
    nbColors=size(colorMap,1);
    if nargin<3 || isempty(range),
        range=1:nbColors;
    end
    if numel(range)==1,
        % only end is specified
        range=[0 range];
    end
    if numel(range)==2,
        % only start and end is given
        range=range(1)+diff(range)*(([1:nbColors]-1)./nbColors);
    end
    if nargin<4 || isempty(defaultValue),
        defaultValue=NaN;
    end
    
    % convert to double to avoid clipping with uints
    imgRGB=double(imgRGB);
    colorMap=double(colorMap(:,:));
    
    sqdDistances=ones(imgSize(1:2))*max(imgRGB(:))*(1+imgSize(3)); % larger than any other possible value
    values=ones(imgSize(1:2));
    for colIdx=1:nbColors,
        color=colorMap(colIdx,:);
        newSqdDistance=sum(abs(imgRGB-repmat(reshape(color,[1 1 3]),imgSize(1:2))).^2,3);
        betterPixels=newSqdDistance<sqdDistances;
        sqdDistances(betterPixels)=newSqdDistance(betterPixels);
        values(betterPixels)=colIdx;
    end
    % Remove outliers
    sqdDistanceInMap=sum(abs(diff(colorMap)).^2,2);
    outliers=sqdDistances>2*max(sqdDistanceInMap);
    values(outliers)=defaultValue;
    % convert to value
    values(~outliers)=range(values(~outliers));
end