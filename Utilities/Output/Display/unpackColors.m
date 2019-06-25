% img=unpackColors(packedColors)
%
% Takes an 2D matrix of uint32 values packed as
% Blue, Green, Red, and Alpha with each uint8 bit values.
% Outputs a 3D matrix with 4 channels: Red, Green, Blue, and Alpha
% with double precision values between 0 and 1.
%
function img=unpackColors(packedColors)
    imgSize=size(packedColors);
    img=typecast(uint32(packedColors(:)),'uint8'); % B G R A
    img=permute(reshape(img,[4 imgSize]),[2 3 1]);
    img(:,:,1:3)=img(:,:,3:-1:1); % change to R G B A
    img=double(img)./255;    
end