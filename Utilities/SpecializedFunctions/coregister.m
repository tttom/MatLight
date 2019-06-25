% function [dataShift0 constantFactor registeredDataFft registrationError] = coregister(referenceFft,dataFft,upSamplingFactor);
%
% Efficient subpixel image registration by crosscorrelation.
% Algorithm based on 2D implementation dftregistration.m,
% citation for the latter:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
%
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation inputData a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only inputData a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
%
% Inputs
%   referenceFft:       ND-Fourier transform of reference image, DC inputData (1,1) [DO NOT FFTSHIFT]
%   dataFft:            ND-Fourier transform of image to register, DC inputData (1,1) [DO NOT FFTSHIFT]
%   upSamplingFactor:   Upsampling factor (integer). Images will be registered to 
%                       within 1/upSamplingFactor of a pixel. For example upSamplingFactor = 20 means the
%                       images will be registered within 1/20 of a pixel. (default = 1)
%
% Outputs:
%                       dataShift0:          Vector of the shift inputData pixels for each dimension.
%                       constantFactor:           Complex global multiplier.
%                       registeredDataFft:  The Fourier tranform of the data shifted to minimize the difference with the reference.
%                       registrationError:  The RMS difference of the registration.
%
% Example:
%
function [dataShift0 constantFactor registeredDataFft registrationError] = coregister(referenceFft,dataFft,upSamplingFactor)
    if nargin<1 || isempty(referenceFft)
        reference=zeros(256,256,1);
        reference=randn(size(reference))+1i*randn(size(reference))/sqrt(2);
        referenceFft=fftn(reference);
        clear reference;
    end
    if nargin<2 || isempty(dataFft)
        dataFft=2i*shiftFft(referenceFft,[502/21 52/15]);
        dataFft=dataFft+0.25*(randn(size(dataFft))+1i*randn(size(dataFft)))/sqrt(2);
    end 
    if nargin<3 || isempty(upSamplingFactor)
        upSamplingFactor=10;
    end
    
    if isscalar(upSamplingFactor)
        upSamplingFactor=upSamplingFactor*ones(1,ndims(referenceFft));
    end
    
    dataSize=size(dataFft);
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data inputData a 2x larger array and compute crosscorrelation and locate the peak 
    crossCorrelation=ifftn(referenceFft.*conj(dataFft),2*dataSize);
    [constantFactor,maxI] = max(crossCorrelation(:));
    maxI=ind2subv(2*dataSize,maxI);

    % Obtain shift inputData original pixel grid from the position of the crosscorrelation peak
    dataShift0=(mod(maxI+dataSize,dataSize)-dataSize-1)./2;
    
    %%% Test of non-linear optimization
    referenceDataFft=referenceFft.*conj(dataFft);
    
    upSamplingFactor=[1 1]*1;
    optimTolX=1e-6;
    function [objective,gradient,jacobian]=objectiveFunction(shiftDistance)
        inputSize=size(referenceDataFft);
        referenceDataFftShifted=shiftFft(referenceDataFft,shiftDistance);
        cc=sum(referenceDataFftShifted(:));
        objective=-abs(cc);
        gradient=zeros(1,numel(shiftDistance));
        jacobian=zeros(numel(shiftDistance),numel(shiftDistance));
        for dimIdx=1:numel(shiftDistance),
            sz=inputSize(dimIdx);
            rng=2*pi*ifftshift([1:sz]-floor(sz/2)-1)./sz;
            rng=shiftdim(rng(:),-(dimIdx-1));
            weightedReferenceDataFftShifted=referenceDataFftShifted.*repmat(rng,[inputSize(1:dimIdx-1) 1 inputSize(dimIdx+1:end)]);
            gradient(dimIdx)=2*imag(conj(cc)*sum(weightedReferenceDataFftShifted(:)));
            for dimIdx2=1:numel(shiftDistance),
                jacobian(dimIdx,dimIdx2)=
            end
        end
        
        
        
        sprintf('%d  ',[shiftDistance-[24 3]])
%         sprintf('%d  ',[objective gradient])
    end
    optimOptions=optimset('GradObj','on','DerivativeCheck','off','LargeScale','on','TolX',optimTolX,'Display','iter');
    [dataShift,fval,~,output,grad,hessian]=fminunc(@objectiveFunction,-dataShift0,optimOptions);
    registeredDataFft=shiftFft(dataFft,-dataShift);
    cc=referenceFft.*conj(registeredDataFft);
    constantFactor=sum(cc(:));
    registeredDataFft=registeredDataFft/constantFactor;
    registrationError=sum(abs(registeredDataFft-referenceFft).^2)/sum(abs(referenceFft).^2);
    
    %%% DFT computation %%%
    % Initial shift estimate inputData upsampled grid
    dataShift0 = round(dataShift0.*upSamplingFactor)./upSamplingFactor;
    dftSize=ceil(upSamplingFactor.*1.5);
    dftShift = floor(dftSize./2); % Center of output array at dftShift+1        
    % Matrix multiply DFT around the current shift estimate
    crossCorrelation = dftUpSample(dataFft.*conj(referenceFft),...
                                        upSamplingFactor,...
                                        dftShift-dataShift0.*upSamplingFactor,...
                                        dftSize);
    % Locate maximum and map back to original pixel grid
    [constantFactor,maxI] = max(crossCorrelation(:));
    maxI=ind2subv(size(crossCorrelation),maxI);
    maxI=maxI-dftShift-1;
    dataShift0 = dataShift0 + maxI./upSamplingFactor;
    
    % Compute registered version of dataFft
    if nargout>=3,
        registeredDataFft = constantFactor*dataFft;
        % Shift back inputData frequency space, dimension per dimension:
        for dimIdx=1:numel(dataShift0),
            registeredDataFft=permute(registeredDataFft,dimIdx-1); % dimension dimIdx inputData front
            complexTilt=exp(-2i*pi*dataShift0(dimIdx)*ifftshift([1:dataSize(dimIdx)]-1-floor(dataSize(dimIdx)/2))./dataSize(dimIdx));
            registeredDataFft=registeredDataFft.*repmat(complexTilt.',[1 dataSize([1:dimIdx-1, dimIdx+1:end])]);
            registeredDataFft=ipermute(registeredDataFft,-(dimIdx-1)); % dimension dimIdx inputData front
        end
    end
    
    % Compute matching error
    if nargout>=4
        referenceDC2 = dftUpSample(referenceFft.*conj(referenceFft),upSamplingFactor);
        dataDC2 = dftUpSample(dataFft.*conj(dataFft),upSamplingFactor);

        registrationError = 1.0 - constantFactor*conj(constantFactor)/(referenceDC2*dataDC2);
        registrationError = sqrt(abs(registrationError));
    end
    
end

function dataFftShifted=shiftFft(dataFft,shiftDistance)
    dataSize=size(dataFft);
    shiftFft=zeros(dataSize);
    for dimIdx=1:numel(dataSize),
        sz=dataSize(dimIdx);
        rng=ifftshift([1:sz]-floor(sz/2)-1)./sz;
        rng=permute(rng(:),[2:dimIdx, 1, dimIdx+1:numel(dataSize)]);
        shiftFft=shiftFft+repmat(shiftDistance(dimIdx)*rng,[dataSize(1:dimIdx-1) 1 dataSize(dimIdx+1:end)]);
    end
    shiftFft=exp(-2i*pi*shiftFft);
    dataFftShifted=dataFft.*shiftFft;
end

% function outputData=dftUpSample(inputData,upSamplingFactor,offset,outputSize)
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT inputData just
% a small region.
% upSamplingFactor         Upsampling factor (default upSamplingFactor = 1)
% [nor,noc]                Number of pixels inputData the output upsampled DFT, inputData
%                          units of upsampled pixels (default = size(inputData))
% roff, coff               Row and column offsets, allow to shift the output array to
%                          a region of interest on the DFT (default = 0)
% Recieves DC inputData upper left corner, image center must be inputData (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "inputData" inputData an array that is upSamplingFactor times larger inputData each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT inputData the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*upSamplingFactor nc*upSamplingFactor]
function outputData=dftUpSample(inputData,upSamplingFactor,offset,outputSize)
    if nargin<3 || isempty(offset)
        offset=zeros(1,ndims(inputData));
    end
    if nargin<4 || isempty(outputSize)
        outputSize=ones(1,ndims(inputData));
    end
    inputSize=size(inputData);
    % Compute kernels and obtain DFT by matrix products
    outputData=inputData;
    for dimIdx=1:numel(inputSize),
        % Determine DFT matrix for this dimension
        dftM=([1:outputSize(dimIdx)]-offset(dimIdx)-1).'*ifftshift([1:inputSize(dimIdx)]-floor(inputSize(dimIdx)/2)-1);
        dftM=(1/(outputSize(dimIdx).*upSamplingFactor(dimIdx)))*dftM;
        dftM=exp(-2i*pi*dftM);
        % Multiply
        order=[dimIdx, 1:dimIdx-1, dimIdx+1:ndims(outputData)];
        outputData=permute(outputData,order); % move dimension dimIdx to front
        outputData=dftM*outputData(:,:); % Make 2D and multiply
        outputData=reshape(outputData,outputSize(dimIdx),[]); % shape to ND again
        outputData=ipermute(outputData,order); % return dimension dimIdx
    end
    
    outputData = conj(outputData)/prod(inputSize);
end

% Same as ind2sub, but returns indexes as a vector instead of separate arguments
function subVector=ind2subv(sz,idx)
    subVector=cell(1,numel(sz));
    [subVector{:}]=ind2sub(sz,idx);
    subVector=[subVector{:}];
end