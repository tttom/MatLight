% RGB = cmyk2rgb(CMYK)
%
% Converts a CMYK image to an RGB image.
% The color is represented in the final dimension of the array.
%
% see also: CMYK = rgb2cmyk(RGB)
%
function RGB = cmyk2rgb(CMYK)
  dataSize = size(CMYK);
  while dataSize(end) == 1
    dataSize = dataSize(1:end-1);
  end
  nbDims = numel(dataSize);
  if dataSize(nbDims) ~= 4
    error('Input to cmyk2rgb must be an array with a final dimension of size 4, not %d!', dataSize(nbDims));
  end
  CMY = reshape(CMYK, [prod(dataSize(1:nbDims-1)) dataSize(nbDims)]);
  K = reshape(CMY(:, 4), [dataSize(1:nbDims-1) 1 1]);
  CMY = reshape(CMY(:, 1:3), [dataSize(1:nbDims-1) 3 1]);
  RGB = bsxfun(@times, 1 - CMY, 1 - K);
end