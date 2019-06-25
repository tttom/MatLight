% [convolved varargout] = spatiallyVariantConvolution(data,kernels,complexTranspose)
%
% Convolve image with bi-linearly interpolated PSF
% The PSFs should be centered at 1+([1:nbKernels]-0.5)*dataSize./nbKernels
%
% Input:
%    data: an n-d array to be convolved
%    kernels: a 2n-d array of kernels to be used, the dimension has to be
%             double that of data, but the size can be anything smaller than the
%             size of data. The first n dimensions represent the pixels of an
%             individual kernel. The final n dimensions represent the tile-indexes
%             with a kernel per tile. The tiles are uniformely spaced to
%             cover the input data size.
%    complexTranspose: a flag indicating if we want the complex transpose
%    of this operation.
%    
%
% Output:
%    convolved: an array the same size as data, the spatially variant convolution of data
%    varargout: cellPosRanges: a cell array indicating the floating point positions of the kernels
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Leave this not in if reusing this code or parts of it. Thanks.
%
%
% Example:
%       nbKernels=[4 6];
%       kernelSize=[1 1]*32;
%       kernels=zeros([kernelSize nbKernels]);
%       [X Y]=ndgrid([1:kernelSize(1)]-floor(kernelSize(1)/2)-1,[1:kernelSize(2)]-floor(kernelSize(2)/2)-1);
%       R2=X.^2+Y.^2;
%       sigma=.2;
%       for kernelIdx=1:prod(nbKernels),
%          kernel=exp(-R2./(2*sigma^2));
%          kernel=kernel./sum(kernel(:));
%          kernels(:,:,kernelIdx)=kernel;
%          sigma=sigma+.2;
%       end
%       [convolved XPos YPos]=spatiallyVariantConvolution(getTestImage('boats'),kernels);
%       imagesc(convolved);
%
function [convolved, varargout] = spatiallyVariantConvolution(data, kernels, complexTranspose)
    if nargin<1
%         data=ones(128,64);
        data=getTestImage('boats');
    end
    if nargin<2
        nbKernels=[4 6];
        kernelSize=[1 1]*32;
        kernels=zeros([kernelSize nbKernels]);
        [X, Y]=ndgrid([1:kernelSize(1)]-floor(kernelSize(1)/2)-1,[1:kernelSize(2)]-floor(kernelSize(2)/2)-1);
        R2=X.^2+Y.^2;
        sigma=.2;
        for kernelIdx=1:prod(nbKernels)
            kernel=exp(-R2./(2*sigma^2));
            kernel=kernel./sum(kernel(:));
            kernels(:,:,kernelIdx)=kernel;
            sigma=sigma+.2;
        end
        clear nbKernels kernel kernelSize X Y R2 sigma kernelIdx;
    end
    if nargin<3 || isempty(complexTranspose)
        complexTranspose=false;
    else
        complexTranspose=true;
    end
    dataSize=size(data);
    kernelsSize=size(kernels);
    nbDims=max(numel(dataSize),ceil(numel(kernelsSize)/2));
    kernelSize=kernelsSize(1:nbDims);
    nbKernels=kernelsSize(nbDims+1:end);
    
    dataSize(end+1:nbDims)=1;
    kernelSize(end+1:nbDims)=1;
    nbKernels(end+1:nbDims)=1;
    
    cellSize=dataSize./nbKernels; % floating point number
    for dimIdx=1:nbDims
        cellPosRanges{dimIdx}=1+cellSize(dimIdx)*([1:nbKernels(dimIdx)]-0.5);
    end
    
    % Loop over all cells, convolving a patch the size of convolutionSize
    convolved=zeros(dataSize);
    for cellIdx=1:prod(nbKernels)
        cellIndexes=ind2subs(nbKernels,cellIdx);
        % Determine which cell to handle now
        for dimIdx=1:nbDims
            rngLims=cellPosRanges{dimIdx}(cellIndexes(dimIdx))+[-1 1]*cellSize(dimIdx)...
                +[-floor(kernelSize(dimIdx)/2) floor((kernelSize(dimIdx)-1)/2)];
            ranges{dimIdx}=[max(1,floor(rngLims(1))):min(dataSize(dimIdx),ceil(rngLims(2)))];
            % Determine the linear weighting factors
            relPositions{dimIdx}=(ranges{dimIdx}-cellPosRanges{dimIdx}(cellIndexes(dimIdx)))./cellSize(dimIdx);
        end
        % set the first and last weights to 1 as if the same kernel is
        % used beyond the edges of the volume
        for dimIdx=1:numel(cellIndexes)
            linearWeights{dimIdx}=max(0,1-abs(relPositions{dimIdx}));
            if cellIndexes(dimIdx)==1
                linearWeights{dimIdx}(relPositions{dimIdx}<0)=1;
            end
            if cellIndexes(dimIdx)==nbKernels(dimIdx)
                linearWeights{dimIdx}(relPositions{dimIdx}>0)=1;
            end
        end
        clear relPositions
        % Do the actual convolution on the cell
        kernel=kernels([1:prod(kernelSize)]+prod(kernelSize)*(cellIdx-1));
        kernel=reshape(kernel,kernelSize);
        if ~complexTranspose
            convolvedSection=convolveSection(data,kernel,ranges,linearWeights);
        else
            convolvedSection=convolveSection(data,kernel,ranges);
        end
        clear kernel;
        if complexTranspose
            convolvedSection=weight(convolvedSection,linearWeights{:});
        end
            
        % Add the cell to the combined results
        convolved(ranges{:})=convolved(ranges{:})+convolvedSection;
    end
    
    
    %
    % Output
    %
    if nargout<1
        [XCellPos, YCellPos]=ndgrid(cellPosRanges{:});
        figure();
        axs(1)=subplot(1,3,1);
        showImage(data,-1); axis image; title('input');
        hold on; scatter(YCellPos(:),XCellPos(:),'filled'); colorbar();
        axs(2)=subplot(1,3,2);
        showImage(convolved,-1); axis image; title('convolved');
        hold on; scatter(YCellPos(:),XCellPos(:),'filled'); colorbar();
        axs(3)=subplot(1,3,3);
        showImage(abs(convolved-data),-1); axis image; title('convolved-data');
        hold on; scatter(YCellPos(:),XCellPos(:),'filled'); colorbar();
        linkaxes(axs);
        
        clear convolved;
    end
    varargout = cellPosRanges;
end

%
% xWeighted=weight(x,varargin)
%
% Returns an array of xWeighted of the same size as the input x and with
% each element weighted by the product of the elements in varargin. 
% varargin must be a set of vectors with a number of elemets corresponding
% to the dimensions of x. 
%
function x = weight(x,varargin)
    inputSize = size(x);
    nbDims = numel(varargin);
    for dimIdx = 1:nbDims
        wts = varargin{dimIdx};
        x = x(:,:);
        x = x.*repmat(wts(:),[1 size(x,2)]);
        x = reshape(x,inputSize);
        x = permute(x,[2:nbDims,1]);
        inputSize = inputSize([2:nbDims,1]);
    end
end
            
            
%
% convolved = convolveSection(data,kernel,varargin)
%
% Convolves the data with the kernel and returns a sub cube only. The data
% array is extended outward if required.
% 
% Inputs:
%     data: the data to be convolved
%     kernel: the convolution kernel
%     ranges: N vectors with integers, where N is equal to the number of 
%         dimensions of data and kernel.
%     linearWeights: if specified, the transposed kernels are used.
%
% Output:
%     convolved: an array of size indicated by ranges
%
function convolved = convolveSection(data,kernel,ranges,linearWeights)
    dataSize=size(data);
    kernelSize=size(kernel);
    if nargin<4
        linearWeights=[];
    end
    complexTranspose=isempty(linearWeights);
    
    startIndexes=cellfun(@(c) c(1),ranges,'UniformOutput',true);
    sectionSize=cellfun(@(c) numel(c),ranges,'UniformOutput',true);
    nbDims=max(numel(sectionSize),numel(kernelSize));
    % Add singleton dimensions to the sizes
    sectionSize(end+1:nbDims)=1;
    kernelSize(end+1:nbDims)=1;
    dataSize(end+1:nbDims)=1;
    
    % Prepare the range to select that data before and after convolution
    convolutionSize=sectionSize+kernelSize-1;
    extensionLength=convolutionSize-sectionSize;
    extensionLengthBegin=floor((1+extensionLength)./2);
    extensionLengthEnd=extensionLength-extensionLengthBegin;
    extendedRanges={};
    extendedWeights={};
    validRanges={};
    for dimIdx=1:nbDims
        % Extend the ranges both sides, to length convolutionSize(dimIdx)
        extendedRanges{dimIdx}=[1:convolutionSize(dimIdx)]-1+startIndexes(dimIdx)-extensionLengthBegin(dimIdx);
        if ~complexTranspose
            % extend the weights to the same length
            extendedWeights{dimIdx}=linearWeights{dimIdx}([ones(1,extensionLengthBegin(dimIdx)),1:end,end*ones(1,extensionLengthEnd(dimIdx))]);
        end
        % clip indexes so to extend the input data if needed
        extendedRanges{dimIdx}=max(1,min(dataSize(dimIdx),extendedRanges{dimIdx}));
        % also keep track of which data is valid
        validRanges{dimIdx}=extensionLengthBegin(dimIdx)+[1:sectionSize(dimIdx)];
    end
    
    % Extend and convolve
    if ~complexTranspose
        selectedData=weight(data(extendedRanges{:}),extendedWeights{:});
    else
        extendedData=data(extendedRanges{:});
        % replace extension by zeros
        selectedData=zeros(size(extendedData));
        selectedData(validRanges{:})=extendedData(validRanges{:});
        clear extendedData;
    end
    convolved=convolveFt(selectedData,kernel,complexTranspose);
    clear selectedData;
    % Crop to original size
    if ~complexTranspose
        convolved=convolved(validRanges{:});
    else
        % integrate the pixels in the extended region
        dataSize=size(convolved);
        for dimIdx=1:nbDims
            sumBegin=sum(convolved(1:(validRanges{dimIdx}(1)-1),:),1);
            sumEnd=sum(convolved((validRanges{dimIdx}(end)+1):end,:),1);
            convolved=convolved(validRanges{dimIdx},:);
            convolved(1,:)=convolved(1,:)+sumBegin;
            convolved(end,:)=convolved(end,:)+sumEnd;
            dataSize(1)=sectionSize(dimIdx); % cropped now
            convolved=reshape(convolved,dataSize);
            convolved=permute(convolved,[2:nbDims,1]);
            dataSize=dataSize([2:nbDims,1]);
        end
    end
end

%
% convolved=convolveFt(data,kernel)
%
% Convolves the data, wrapping it, not extending the data.
%
% Inputs:
%     data: the data to be convolved
%     kernel: the convolution kernel
%     complexTranspose: if true, the kernel is flipped in all dimensions,
%                       default: false
%
% Output: convolved, a matrix of the same size as data
%
% data and kernel have to have the same number of dimensions but not the
% same dimensions. The size of data should be at least as large as the
% kernel
%
function convolved=convolveFt(data,kernel,complexTranspose)
    if nargin<3
        complexTranspose=false;
    end
    dataSize=size(data);
    kernelSize=size(kernel);
    nbDims=max(numel(dataSize),numel(kernelSize));
    % Add singleton dimensions to the sizes
    dataSize(end+1:nbDims)=1;
    kernelSize(end+1:nbDims)=1;
    
    convolutionSize=max(dataSize,kernelSize);
    % zero pad the kernel
    if any(kernelSize<convolutionSize)
        % zero pad and center on first pixel
        convolutionSizeCell=num2cell(convolutionSize);
        kernel(convolutionSizeCell{:})=0;
    end
    % Center on the top right element as ifftshift
    kernel=circshift(kernel,-floor(kernelSize./2));
    % Convolve
    transferFunction=fftn(kernel);
    if complexTranspose
        transferFunction=conj(transferFunction);
    end
    if ~isreal(data) || ~isreal(kernel)
        convolved=ifftn(fftn(data).*transferFunction);
    else
        % faster and avoids rounding errors
        convolved=ifftn(fftn(data).*transferFunction,'symmetric');
    end
end
