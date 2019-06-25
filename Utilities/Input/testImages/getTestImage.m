%
% testImage=getTestImage(fileName)
%
% Reads images from the testImages sub-folder placed next to this function.
% The returned matrix of doubles is scaled so the intensity values are in [0,1].
%
% Example usage:
%   testScene=getTestImage('lena');
%   testScene=getTestImage('boats');
%   testScene=getTestImage('lena512bw.png');
%   testScene=getTestImage('boat.512.tiff');
%   testScene=getTestImage('spoke');%built-in
%
function testImage=getTestImage(fileName)
    switch(fileName)
        case 'spoke'
            testImage = getSpokeTarget([1 1]*512,40,.8);
        case {'boat', 'boats'}
            testImage = getTestImage('boat.512.tiff');
        case 'usaf'
            testImage = getTestImage('usaf1951_750x750.png');
      otherwise
            filePath = fileparts(mfilename('fullpath'));
            fullFile = fullfile(filePath, fileName);
            try
                img = imread(fullFile);
            catch Exc
                img = imread(fileName);
                fullFile = fileName;
            end
            info = imfinfo(fullFile);
            bitsPerColorChannel=info.BitDepth;
            if strcmpi(info.ColorType,'truecolor')
                bitsPerColorChannel=bitsPerColorChannel/3;
            end
            testImage=double(img)/(2^bitsPerColorChannel);
    end    
end