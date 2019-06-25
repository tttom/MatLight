% [dataCube maxPotentialValue]=readDataCubeFromImageFolder(folderName,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     fileName: string representing the folder containing a series of image files numbered starting from 0.png
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%
% Outputs:
%     dataCube: the 3D matrix of values
%     maxPotentialValue: the maximum value that could have been stored in this file.
%
% Returns a 3D matrix in single precision and normalized to 1;
%
function [dataCube, maxPotentialValue]=readDataCubeFromImageFolder(folderName,projectionDimensionOrSubCube,frameIndexes,normalize)
    if nargin<1 || isempty(folderName)
        folderName='.';
    end
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    if (nargin<3)
        frameIndexes=-1;
    end
    frameIndexes=sort(frameIndexes);
    if (nargin<4 || isempty(normalize))
        normalize=true;
    end
    
    % Check if the specified folder is actually a file
    if ~isempty(regexpi(folderName,'\.(png|tiff?|jpe?g)$'))
        [folderName, fileName, extension]=fileparts(folderName);
        imageFileRegexp=strcat('^',regexprep(fileName,'\*','(?<nb>\\d+)','ignorecase'),'\',extension);
    else
        imageFileRegexp='(?<nb>\d+)\.(png|tiff?|jpe?g)'; %anything image with a file name ending in a number
    end
    
    allImageFiles=dir(fullfile(folderName,'*'));
    allImageFiles={allImageFiles.name};
    imageFiles={};
    for fileIdx=1:numel(allImageFiles)
        if ~isempty(regexpi(allImageFiles{fileIdx},imageFileRegexp))
            imageFiles{end+1}=allImageFiles{fileIdx};
        end
    end
    clear allImageFiles;
    
    % Find out how many frames we have
    nbFrames=numel(imageFiles);
    
    % sort ascending by the numbers in the file name, left most is most significant
    indices=zeros(nbFrames,1); % row: current frame index, col: file name numbers
    for fileIdx=1:nbFrames
        tokens=regexpi(imageFiles{fileIdx},imageFileRegexp,'tokens');
        nbIdx=1;
        for tokenIdx=1:numel(tokens)
            for tokenIdxIdx=1:numel(tokens{tokenIdx})
                numberString=regexp(tokens{tokenIdx}{tokenIdxIdx},'\d+','match');
                if ~isempty(numberString)
                    indices(fileIdx,nbIdx)=str2double(numberString);
                    nbIdx=nbIdx+1;
                end
            end
        end
    end
    [~,sortedI]=sortrows(indices);
    imageFiles=imageFiles(sortedI);
    clear indices
    
    getFrameFileName=@(frameIdx) fullfile(folderName,imageFiles{frameIdx});
    function img = getFrame(frameIdx)
        fileName = getFrameFileName(frameIdx);
        try
            img = mean(single(imread(fileName)),3);
        catch Exc
            logMessage('Failed to read file %s!',fileName);
            throw(Exc);
        end
    end
    
    if (numel(frameIndexes)==1 && frameIndexes==-1)
        frameIndexes=[1:nbFrames];
    else
        frameIndexes=intersect([1:nbFrames],frameIndexes);
    end

    subCubeSliced = max(size(projectionDimensionOrSubCube))>1;
    if subCubeSliced
        if size(projectionDimensionOrSubCube,1)>=3 % also sliced in the third dimension
            frameIndexes = frameIndexes(frameIndexes>=projectionDimensionOrSubCube(3,1) ...
                                        & frameIndexes<=projectionDimensionOrSubCube(3,2));
        end
    end
    
    if ~isempty(frameIndexes)
        firstFrame = getFrame(frameIndexes(1));
        imgSize = size(firstFrame);
        if subCubeSliced
            imgSize = min(imgSize,1+diff(projectionDimensionOrSubCube(1:2,:).'));
        end

        imgInfo = imfinfo(getFrameFileName(frameIndexes(1)));
        maxPotentialValue = 2^imgInfo.BitDepth-1;
        
        dataCube=zeros([imgSize numel(frameIndexes)],'single');
        for frameIndexIdx = 1:numel(frameIndexes)
            frameIndex = frameIndexes(frameIndexIdx);
            if frameIndexIdx>1
                img=getFrame(frameIndex);
            else
                img=firstFrame; % don't read it twice
            end
            % (project and) store
            if isempty(projectionDimensionOrSubCube)
                %Complete data cube
                dataCube(:,:,frameIndexIdx)=img;
            else
                %Project or crop the data cube
                if (max(size(projectionDimensionOrSubCube))==1)
                    %Project the full cube along one dimension specified by projectionDimensionOrSubCube
                    if (any(projectionDimensionOrSubCube==1))
                        img=max(img,[],1);
                    end
                    if (any(projectionDimensionOrSubCube==2))
                        img=max(img,[],2);
                    end
                    if (~any(projectionDimensionOrSubCube==3))
                        dataCube(:,:,frameIndexIdx)=img;
                    else
                        dataCube(:,:,1)=dataCube(:,:,1)+img;
                    end
                else
                    %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
                     img = img(projectionDimensionOrSubCube(1,1):projectionDimensionOrSubCube(1,2),projectionDimensionOrSubCube(2,1):projectionDimensionOrSubCube(2,2));
                     dataCube(:,:,frameIndexIdx) = img;
                end
            end
        end
    
        if (normalize)
            dataCube=dataCube./maxPotentialValue;
        end
    else
        logMessage('No data found in %s!', folderName);
        dataCube=[];
    end
end