% [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateBlurAndSNR(img,xRange,yRange,nbCells,minSubImageSize)
%
% Determines the image blur width and signal to noise assuming the input image img is blurred by a Gaussian kernel and
% independent noise such as Guassian or Poisson has been added during recording.
% To determine the noise it is assumed that the blur is circular and that
% the detector array oversamples.
% The signal is modeled as inversely proportional to the spatial frequency,
% blurred with a Gaussian and noise added.
% The fitting error is measuredusing the l2-norm and the error is weighted
% with the spatial frequency.
%
% x/yRange: optional: the coordinates of the rows and columns.
% nbCells: optional: a vector of how the image will be subdivided, default
%                    no subdivisions
% minSubImageSize: optional: the minimum size of the sub images to be
% used, default: ceil(size(img)./nbCells)
%
% Output: abs(fft2(img)) ~ (amplitude/f)exp(-f^2/(2*sigmaFt^2)) + noiseLevelFt
%    where sigmaFt=1./(sigma*2*pi) and noiseLevelFt=noiseLevel*sqrt(prod(imgSize(1:2)))
%    x/yfRange: The spatial frequency ranges corresponding to x/yRange
%
% Example:
%     groundTruth=getTestImage('boats');
%     imgSize=size(groundTruth);
%     xRange=[1:imgSize(1)]; yRange=[1:imgSize(2)];
%     sigma=5;
%     noiseLevel=0.05;
%     [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange);
%     [XF,YF]=ndgrid(xfRange,yfRange);
%     RF2=XF.^2+YF.^2;
%     clear XF YF;
%     sigmaFt=1./(sigma*2*pi);
%     GFt=exp(-RF2./(2*sigmaFt^2));
%     img=ifft2(fft2(groundTruth).*repmat(GFt,[1 1 imgSize(3:end)]),'symmetric');
%     img=img+noiseLevel*randn(size(img)); % additive
%     clear sigma noiseLevel xfRange yfRange RF2 sigmaFt GFt
%     [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateBlurAndSNR(img);
%     [XF,YF]=ndgrid(xfRange,yfRange);
%     RF=sqrt(XF.^2+YF.^2);
%     clear XF YF
%     otf=exp(-RF.^2./(2*sigmaFt^2));
% %     amplitude=amplitude/2; % be more conservative
%     NSR=RF.*(noiseLevelFt/amplitude);
%     filter=conj(otf)./(abs(otf).^2+NSR.^2);
%     restoredImg=ifft2(fft2(img).*ifftshift(filter),'symmetric');
%     close all; figure();
%     axs(1)=subplot(1,2,1);
%     showImage(img);
%     axs(2)=subplot(1,2,2);
%     showImage(restoredImg);
%     linkaxes(axs);
%     
function [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateBlurAndSNR(img,xRange,yRange,nbCells,minSubImageSize)
    if nargin<1 || isempty(img),
        randn('seed',0);
        img=getTestImage('boats');
%         img=getTestImage('C:\Users\tvettenb\Downloads\test.jpg');
        
        imgSize=size(img); imgSize(end+1)=1;
        xRange=[1:imgSize(1)]; yRange=[1:imgSize(2)];
        GFt=gaussian2DFt(xRange,yRange,2);
        img=ifft2(fft2(img).*repmat(GFt,[1 1 imgSize(3:end)]),'symmetric');
        noiseLevel=0.01;
%         img=img.*(1+noiseLevel*randn(size(img))); % Poisson like
        img=img+noiseLevel*randn(size(img)); % additive
        clear xRange yRange imgSize noiseLevel GFt
    end
    imgSize=size(img); imgSize(end+1:3)=1;
    if nargin<2 || isempty(xRange),
        xRange=[1:imgSize(1)];
    end
    if nargin<3 || isempty(yRange),
        yRange=[1:imgSize(2)];
    end
    if nargin<4 || isempty(nbCells),
        nbCells=imgSize*0+1;
    end
    if nargin<5 || isempty(minSubImageSize),
        minSubImageSize=[];
    end
    minSubImageSize(end+1:3)=1;
    
    samplePitch=[diff(xRange(1:2)) diff(yRange(1:2))];
    
    % Slice the input image and process cell by cell
    cellSize=imgSize./nbCells;
    if isempty(minSubImageSize),
        subImageSize=ceil(cellSize);
    else
        subImageSize=max(ceil(cellSize),minSubImageSize);
    end
    xCellRange=1+floor(cellSize(1)*([1:nbCells(1)]-0.5));
    yCellRange=1+floor(cellSize(2)*([1:nbCells(2)]-0.5));
    zCellRange=1+floor(cellSize(3)*([1:nbCells(3)]-0.5));
    [XCell,YCell,ZCell]=ndgrid(xCellRange,yCellRange,zCellRange);
    xRangeSel=samplePitch(1)*[1:subImageSize(1)]; yRangeSel=samplePitch(2)*[1:subImageSize(2)];
    sigma=zeros(nbCells); amplitude=zeros(nbCells); noiseLevel=zeros(nbCells); sigmaFt=zeros(nbCells); noiseLevelFt=zeros(nbCells); relFitError=zeros(nbCells);
    for idx=1:prod(nbCells),
        imgSel=circshift(img,1+floor(subImageSize./2)-[XCell(idx) YCell(idx) ZCell(idx)]);
        imgSel=imgSel(1:subImageSize(1),1:subImageSize(2),1:subImageSize(3));
        [sigma(idx) amplitude(idx) noiseLevel(idx) xfRange yfRange sigmaFt(idx) noiseLevelFt(idx) relFitError(idx)]=estimateBlurAndSNRSingleImage(imgSel,xRangeSel,yRangeSel);
    end
end
function [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateBlurAndSNRSingleImage(img,xRange,yRange)
    samplePitch=[diff(xRange(1:2)) diff(yRange(1:2))];
    imgSize=size(img); imgSize(end+1:3)=1;
    
    % Calculate the image spectrum
    imgSpectrum=fftshift(sqrt(mean(abs(fft2(img)).^2,3)));
    clear img;
        
    % Model the image content in the spatial frequency domain 
    [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange,'centered');
    [XF,YF]=ndgrid(xfRange,yfRange);
    [~,RF]=cart2pol(XF,YF);
    reliableFreqs=(abs(XF)>0 & abs(YF)>0); % due to wrapping discontinuities, artificial horizontal and vertical freqs may bias the result
%     reliableFreqs=(abs(XF)>diff(xfRange(1:2)) & abs(YF)>diff(xfRange(1:2))); % do a bit more just in case
    clear XF YF;
    unitySpectrum=calcNaturalSpectrum(xfRange,yfRange); % must be centered
    
    %
    % Estimate noise level
    %
    % Assume at least Nyquist sampling in both dimensions
    % Use only the highest spatial frequencies
    opticalCutOffFreq=min(0.5./samplePitch);
    noiseFreqs=(RF>opticalCutOffFreq) & reliableFreqs;
    noiseLevelFt=sqrt(sum(imgSpectrum(noiseFreqs).^2)/sum(noiseFreqs(:)));
    noiseLevel=noiseLevelFt/sqrt(prod(imgSize(1:2)));
    
    l2Norm=@(x) sqrt(sum(abs(x(:)).^2));
    calcError=@(a) l2Norm((a-imgSpectrum).*RF.*reliableFreqs);
    
    testSpectrum=@(p) sqrt((p(1)*unitySpectrum.*fftshift(gaussian2DFt(xRange,yRange,p(2)))).^2+noiseLevelFt^2);
%     testSpectrum=@(p) max(p(1)*unitySpectrum.*fftshift(gaussian2DFt(xRange,yRange,p(2))),noiseLevelFt);
    
    p0=[1.0 1.0]; % [sqrt(amplitude) sqrt(sigma)]
    [pOpt fitError]=fminsearch(@(p) calcError(testSpectrum(abs(p).^2)),p0,optimset('Display','none','TolX',1e-1));
    amplitude=abs(pOpt(1))^2;
    sigma=abs(pOpt(2)).^2;
    relFitError=fitError/l2Norm(imgSpectrum.*RF);
    
    sigmaFt=1./(sigma*2*pi);
        
    if nargout<1,
        logMessage('The img with pixel pitch (%0.3d,%0.3d) is blurred with a standard deviation %0.3f: ~exp(-p^2/(2*%0.3f^2)), and has %0.3f%% noise.',[samplePitch sigma sigma 100*noiseLevel]);
        logMessage('The magnitude of its Fourier transform abs(fft2(img)) ~ (%0.3f/f)exp(-f^2/(2*%0.3f^2)) + %0.6f.',[amplitude sigmaFt noiseLevelFt]);
    end
end

function [GFt xfRange yfRange]=gaussian2DFt(xRange,yRange,sigma)
    [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange);
    [XF,YF]=ndgrid(xfRange,yfRange);
    RF2=XF.^2+YF.^2;
    clear XF YF;
    sigmaFt=1./(sigma*2*pi);
    GFt=exp(-RF2./(2*sigmaFt^2));
end
