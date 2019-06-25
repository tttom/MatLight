% writeDataCubeToTiffFile(dataCube,fileName)
%
% Stores a 3D matrix with the frames of a tiff file.
%
% Inputs:
%     dataCube: a 3D matrix of real numbers between o and 1.
%     fileName: string representing the tiff file to read as a data cube
%
function writeDataCubeToTiffFile(dataCube,fileName)
    nbFrames=size(dataCube,3);
    
    maxValue=2^16-1;
    
    for frameIdx=1:nbFrames,
        img=dataCube(:,:,frameIdx);
        img=uint16(img*maxValue);
        if (frameIdx>1)
            writeMode='append';
        else
            writeMode='overwrite';
        end
        imwrite(img,fileName,'tiff',...
            'Compression','deflate','WriteMode',writeMode);
    end
end

