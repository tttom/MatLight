% packedColors=packColors(img)
%
% Takes an 3D matrix of 2D images stacked as
% Red, Green, Blue, and optionally, Alpha
% with values between 0 and 1.
% Outputs a 2D matrix with 32 bit values of
% with 8 bits per channel.
%
function packedColors=packColors(img)
    imgSize=size(img);
    img(:,:,3:-1:1)=img(:,:,1:3); % change to B G R (A)
    img(:,:,(end+1):4)=0; % Add alpha channel if not present
    img=uint8(255*img);
    img=permute(img,[3 1 2]);
    packedColors=typecast(img(:),'uint32');
    packedColors=reshape(packedColors,imgSize(1:2));
    %packedColors=double(packedColors);
end