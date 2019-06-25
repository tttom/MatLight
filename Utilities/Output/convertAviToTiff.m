% convertAviToTiff(aviFileName,tiffFileName)
%
% Converts a 3D matrix from avi 16 bit format to 16 bit tiff format (lossless compressed).
%
% Inputs:
%     aviFileName: the avi file name used as input.
%     tiffFileName: (optional) the output file name.
%
function convertAviToTiff(aviFileName,tiffFileName)
    if (~strcmpi(aviFileName(max(1,end-3):end),'.avi'))
        aviFileName=[aviFileName,'.avi'];
    end
    if (nargin<2 || isempty(tiffFileName))
        tiffFileName=[aviFileName(1:end-4), '.tif'];
    end
    
    if (exist('VideoReader','class'))
        vidReader = VideoReader(aviFileName, 'tag', 'vidReader1');

        % Get the data size.
        nbFrames = vidReader.NumberOfFrames;
        
        colorEncodingFor16bits=vidReader.BitsPerPixel==24;
        
        for frameIdx = 1:nbFrames,
            % Read
            vidFrame = read(vidReader,frameIdx);
            if (colorEncodingFor16bits)
                % 16 bit encoded in green and blue channel
                img = uint16(vidFrame(:,:,3))+uint16(vidFrame(:,:,2))*256;
            else
                img = uint8(mean(single(vidFrame),3));
            end
            
            %Write
            if (frameIdx>1)
                writeMode='append';
            else
                writeMode='overwrite';
            end
            
            nbTrials=10;
            while nbTrials>0
                try
                    nbTrials=nbTrials-1;
                    imwrite(img,tiffFileName,'tiff',...
                        'Compression','deflate','WriteMode',writeMode);
                    nbTrials=0;
                catch Exc
                    if (nbTrials>0)
                        logMessage('Attempting to open file for writing.');
                    else
                        logMessage('Error opening file for writing.');
                    end
                    pause(1);
                end
            end
        end
        
    else
        logMessage('VideoReader class not found, Matlab version not appropriate.');
    end
end

