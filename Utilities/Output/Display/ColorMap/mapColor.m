% imgRGB=mapColor(img,colorMap)
%
% Converts a gray scale image with values between zero and 1 to an RGB
% image using a color map.
%
% Input parameters:
%     img: a two dimensional matrix of intensity values, either uint16,
%          uint8, or values between 0 and 1.
%     colorMap: a function handle with one input argument, the intensity,
%               or a matrix with three values per row (Red, Green, or Blue),
%               and one row per uniformely-spaced intensity.
%
% Output: a three-dimensional matrix of floating point values between 0 and
%         one, the third dimension indicating the color channel (Red, Green
%        , or Blue).
function imgRGB=mapColor(img,colorMap)
    if (isa(colorMap,'function_handle'))
        colorMap=colorMap(4096);
    end
    if (isinteger(img))
        switch class(img)
            case 'uint16'
                img=double(img)./(2^16-1);
            otherwise
                %Assume 8 bit
                img=double(img)./(2^8-1);
        end
    end
    
    inputSize = size(img);
    outputSize = [inputSize, size(colorMap,2)];
    
    nbColors = size(colorMap,1);
    img = 1+floor(nbColors*img);
    img = min(max(1,img),nbColors);
    imgRGB = reshape(colorMap(img,:),outputSize);
end