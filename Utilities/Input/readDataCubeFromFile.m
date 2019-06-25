% [dataCube maxValue]=readDataCubeFromFile(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     dataSource: string representing the file to read as a data cube, a matfile object, or a VideoReader object.
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%
% Outputs:
%     dataCube: Returns a 3D matrix with the frames of an avi file in single precision and normalized to 1 unless the normalize argument is false;
%     maxValue: the maximum value that could have been stored in this file.
%
function [dataCube maxValue]=readDataCubeFromFile(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize)
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    if (nargin<3)
        frameIndexes=-1;
    end
    if (nargin<4)
        normalize=true;
    end
    
    showProcessedInMatFile=true;
    
    if (isa(dataSource,'matlab.io.MatFile'))
        if (showProcessedInMatFile && isprop(dataSource,'restoredDataCube'))
            dataCube=dataSource.restoredDataCube(:,:,frameIndexes);
        else
            dataCube=dataSource.recordedImageStack(:,:,frameIndexes);
        end
        maxValue=max(dataCube(:));
        if (normalize)
            dataCube=dataCube./maxValue;
        end
    else
        if (~ischar(dataSource))
            % It must be a VideoReader object
            [dataCube maxValue]=readDataCubeFromAviFile(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize);
        else
            if ~isdir(dataSource)
                fileType=lower(dataSource(max(1,end-3):end));

                switch(fileType)
                    case '.avi'
                        [dataCube maxValue]=readDataCubeFromAviFile(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize);
                    case {'.tif','tiff'}
                        [dataCube maxValue]=readDataCubeFromTiffFile(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize);
                    otherwise
                        logMessage('Error: Unknown file type');
                end
            else
                if exist(fullfile(dataSource,'0.png'),'file')
                    [dataCube maxValue]=readDataCubeFromPngFolder(dataSource,projectionDimensionOrSubCube,frameIndexes,normalize);
                else
                    logMessage('Folder %s does not contain a file with the name ''0.png''.',dataSource);
                end
            end
        end
    end
end

