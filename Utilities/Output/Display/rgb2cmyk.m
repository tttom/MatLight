% CMYK = rgb2cmyk(RGB)
%
% Converts an RGB image to a CMYK image.
% The color is represented in the final dimension of the array.
%
% see also: RGB = cmyk2rgb(CMYK)
%
function CMYK = rgb2cmyk(RGB)
  dataSize = size(RGB);
  while dataSize(end) == 1
    dataSize = dataSize(1:end-1);
  end
  nbDims = numel(dataSize);
  if dataSize(nbDims) ~= 3
    error('Input to rgb2cmyk must be an array with a final dimension of size 3, not %d!', dataSize(nbDims));
  end
  CMY = 1 - RGB;
  K = min(CMY, [], nbDims);
  CMY = bsxfun(@rdivide, bsxfun(@minus, CMY, K), 1 - K);
  CMYK = cat(nbDims, CMY, K);
end