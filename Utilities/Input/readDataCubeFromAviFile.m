% [dataCube maxValue]=readDataCubeFromAviFile(fileName,projectionDimension,dataRanges,normalize,channelsToUse)
%
% Inputs:
%     fileName: string representing the avi file to read as a data cube, or a VideoReader object.
%     projectionDimension: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     dataRanges: optional vector of frame indexes to read (base 1), or a
%                 cell array with index vectors for rows, colums, and frames. An empty
%                 list indicates all available indexes will be read. Default: all read.
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%     channelsToUse: one or more channels Red, Green, Blue, to use. When
%                   this is the empty list, all channels will be used. Default: []
%
% Outputs:
%     dataCube: Returns a 3D matrix with the frames of an avi file in single precision and normalized to 1 unless the normalize argument is false;
%     maxValue: the maximum value that could have been stored in this file.
%
function [dataCube maxValue]=readDataCubeFromAviFile(fileName,projectionDimension,dataRanges,normalize,channelsToUse)
    if (nargin<2)
        projectionDimension=[];
    end
    if (nargin<3)
        dataRanges=[];
    end
    if (nargin<4)
        normalize=true;
    end
    if nargin<5,
        channelsToUse=[];
    end
    
	if (exist('VideoReader','class'))
        if (ischar(fileName))
            vidReader = VideoReader(fileName, 'tag', 'vidReader1');
        else
            vidReader=fileName;
        end

        % Get the data maximum size.
        dataSize = [vidReader.Height, vidReader.Width, vidReader.NumberOfFrames];
        % Convert a frames range to a set of ranges for each dimension
        if ~iscell(dataRanges),
            dataRanges={[],[],dataRanges}; %rowRange, colRange, frameRange
        end
        % keep only existing indexes
        dataRanges=cellfun(@(r) sort(r),dataRanges,'UniformOutput',false);
        for rangeIdx=1:max(numel(dataRanges),3),
            if rangeIdx>numel(dataRanges) || isempty(dataRanges{rangeIdx}),
                dataRanges{rangeIdx}=[1:dataSize(rangeIdx)];
            else
                dataRanges{rangeIdx}=intersect([1:dataSize(rangeIdx)],dataRanges{rangeIdx});
            end
        end
        % Reduce the data maximum size if needed.
        dataSize = cellfun(@(r) numel(r),dataRanges);
    
        colorEncodingFor24bits=vidReader.BitsPerPixel==24;
        
        if isempty(channelsToUse),
            channelsToUse=[1:floor(vidReader.BitsPerPixel/8)];
        end
        nbChannelsToUse=numel(channelsToUse);
        
        % Create a 3D matrix from the video frames.
        dataCube=zeros(dataSize,'single');
        for frameIndexIdx = 1:dataSize(3),
            frameIndex = dataRanges{3}(frameIndexIdx);
            vidFrame = single(read(vidReader,frameIndex));
            vidFrame = vidFrame(dataRanges{1},dataRanges{2},channelsToUse);
            if (colorEncodingFor24bits && nbChannelsToUse>1),
                % 16 bit encoded in green and blue channel
                img=zeros(dataSize(1:2),'single');
                for channelIdx=1:nbChannelsToUse,
                    img = (2^8)*img + vidFrame(:,:,channelIdx);
                end
            else
                img = mean(vidFrame,3);
            end
            if (isempty(projectionDimension))
                %Complete data cube
                dataCube(:,:,frameIndexIdx)=img;
            else
                %Project or crop the data cube
                if (max(size(projectionDimension))==1)
                    %Project the full cube along one dimension specified by projectionDimension
                    if (any(projectionDimension==1))
                        img=max(img,[],1);
                    end
                    if (any(projectionDimension==2))
                        img=max(img,[],2);
                    end
                    if (~any(projectionDimension==3))
                        dataCube(:,:,frameIndexIdx)=img;
                    else
                        dataCube(:,:,1)=dataCube(:,:,1)+img;
                    end
                end
            end
        end
        if (ischar(fileName))
            delete(vidReader);
        end
    else
         warning off;
             mov=aviread(fileName,frameIdx);
         warning on;
         colorEncodingFor24bits=false;
         
         %Convert to matrix
         dataCubeCellArray=struct2cell(mov);
         dataCubeCellArray=dataCubeCellArray(1,1,:);
         dataCube=cell2mat(dataCubeCellArray);
         clear dataCubeCellArray;
	end
    
    if (colorEncodingFor24bits)
        maxValue=2^(8*nbChannelsToUse)-1;
    else
        maxValue=2^8-1;
    end
    if (normalize)
        dataCube=dataCube./maxValue;
    end
end