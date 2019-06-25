% [dataCube maxValue]=readDataCubeFromTiffFile(fileName,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     fileName: string representing the tiff file to read as a data cube
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input file, default: true
%
% Outputs:
%     dataCube: the 3D matrix of values
%     maxValue: the maximum value that could have been stored in this file.
%
% Returns a 3D matrix with the frames of a tiff file in single precision and normalized to 1;
%
function [dataCube, max_value]=readDataCubeFromTiffFile(fileName, projectionDimensionOrSubCube, frameIndexes, normalize)
    if nargin < 2
        projectionDimensionOrSubCube = [];
    end
    if nargin < 3
        frameIndexes = -1;
    end
    frameIndexes = sort(frameIndexes);
    if nargin<4 || isempty(normalize)
        normalize = true;
    end
    
    info = imfinfo(fileName);
    img_shape = [info(1).Height info(1).Width];
    nbFrames = length(info);
%     nbChannels = max([info.SamplesPerPixel]);
    max_value = 2^(info(1).BitDepth)-1;
        
    if length(frameIndexes) == 1 && frameIndexes(1) == -1
        frameIndexes = [1:nbFrames];
    else
        frameIndexes = intersect([1:nbFrames],frameIndexes);
    end
    
    % Calculate the size of the output and allocate space
    output_shape = [img_shape, length(frameIndexes)];
    if ~isempty(projectionDimensionOrSubCube)
      %Project or crop the data cube
      if max(size(projectionDimensionOrSubCube)) == 1
          %Project the full cube along one dimension specified by projectionDimensionOrSubCube
          if any(projectionDimensionOrSubCube == 1)
              output_shape(1) = 1;
          end
          if any(projectionDimensionOrSubCube == 2)
              output_shape(2) = 1;
          end
          if any(projectionDimensionOrSubCube == 3)
              output_shape(3) = 1;
          end
      else
          %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
          output_shape = 1 + diff(projectionDimensionOrSubCube, 1, 2).';
          frameIndexes = frameIndexes(projectionDimensionOrSubCube(3, 1):projectionDimensionOrSubCube(3, 2));
      end
    end
    % Allocate space
    dataCube = zeros(output_shape, 'single');
    % Fill in the data
    for frameIndexIdx = 1:length(frameIndexes)
        frame_index = frameIndexes(frameIndexIdx);
        img = single(imread(fileName, 'Index', frame_index));
        if (normalize)
            img = img ./ single(max_value);
        end
        % (project and) store
        if isempty(projectionDimensionOrSubCube)
            %Complete data cube
            dataCube(:,:,frameIndexIdx) = img;
        else
            %Project or crop the data cube
            if max(size(projectionDimensionOrSubCube)) == 1
                %Project the full cube along one dimension specified by projectionDimensionOrSubCube
                if any(projectionDimensionOrSubCube == 1)
                    img = max(img, [], 1);
                end
                if any(projectionDimensionOrSubCube == 2)
                    img = max(img, [], 2);
                end
                if ~any(projectionDimensionOrSubCube == 3)
                    dataCube(:,:,frameIndexIdx) = img;
                else
                    dataCube(:,:,1) = dataCube(:,:,1) + img;
                end
            else
                %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
                if frame_index >= projectionDimensionOrSubCube(3,1) && frame_index <= projectionDimensionOrSubCube(3,2)
                     img = img(projectionDimensionOrSubCube(1,1):projectionDimensionOrSubCube(1,2),projectionDimensionOrSubCube(2,1):projectionDimensionOrSubCube(2,2));
                     dataCube(:,:,frameIndexIdx) = img;
                end
            end
        end
    end
end

