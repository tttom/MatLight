% [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateGaussianBlurSigma(img,xRange,yRange)
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
%     [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateGaussianBlurSigma(img);
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
function [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateGaussianBlurSigma(img,xRange,yRange)
    if nargin<1 || isempty(img),
%         randn('seed',0);
        img=getTestImage('boats');
%         img=getTestImage('C:\Users\tvettenb\Downloads\test.jpg');
        
        imgSize=size(img); imgSize(end+1)=1;
        xRange=[1:imgSize(1)]; yRange=[1:imgSize(2)];
        GFt=gaussian2DFt(xRange,yRange,3);
        img=ifft2(fft2(img).*repmat(GFt,[1 1 imgSize(3:end)]),'symmetric');
        noiseLevel=0.05;
%         img=img.*(1+noiseLevel*randn(size(img))); % Poisson like
        img=img+noiseLevel*randn(size(img)); % additive
        clear xRange yRange imgSize noiseLevel GFt
    end
    imgSize=size(img);
    if nargin<2 || isempty(xRange),
        xRange=[1:imgSize(1)];
    end
    if nargin<3 || isempty(yRange),
        yRange=[1:imgSize(2)];
    end
    samplePitch=[diff(xRange(1:2)) diff(yRange(1:2))];
    
    % Calculate the image spectrum
%     img=medianFilter(img);
%     img=medianFilter(img);
    imgSpectrum=fftshift(sqrt(mean(abs(fft2(img)).^2,3)));
    clear img;
        
    % Model the image content in the spatial frequency domain 
    [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange,'centered');
    [XF,YF]=ndgrid(xfRange,yfRange);
    [~,RF]=cart2pol(XF,YF);
    clear XF YF;
    unitySpectrum=calcNaturalSpectrum(xfRange,yfRange); % must be centered
    
    %
    % Estimate noise level
    %
    % Assume at least Nyquist sampling in both dimensions
    % Use only the highest spatial frequencies
    opticalCutOffFreq=min(0.5./samplePitch);
    noiseFreqs=RF>opticalCutOffFreq;
    noiseLevelFt=sqrt(sum(imgSpectrum(noiseFreqs).^2)/sum(noiseFreqs(:)));
    noiseLevel=noiseLevelFt/sqrt(prod(imgSize(1:2)));
    
    l2Norm=@(x) sqrt(sum(abs(x(:)).^2));
    calcError=@(a) l2Norm((a-imgSpectrum).*RF);
    
    testSpectrum=@(p) sqrt((p(1)*unitySpectrum.*fftshift(gaussian2DFt(xRange,yRange,p(2)))).^2+noiseLevelFt^2);
%     testSpectrum=@(p) max(p(1)*unitySpectrum.*fftshift(gaussian2DFt(xRange,yRange,p(2))),noiseLevelFt);

    p0=[1.0 1.0]; % sqrt(amplitude) sqrt(sigma)
    [pOpt fitError]=fminsearch(@(p) calcError(testSpectrum(abs(p).^2)),p0,optimset('Display','none','TolX',1e-1));
    amplitude=abs(pOpt(1))^2;
    sigma=abs(pOpt(2)).^2;
    relFitError=fitError/l2Norm(imgSpectrum.*RF);
    
    sigmaFt=1./(sigma*2*pi);
        
    if nargout<1,
        logMessage('The img with pixel pitch (%0.3d,%0.3d) is blurred with a standard deviation %0.3f: ~exp(-p^2/(2*%0.3f^2)), and has %0.3f%% noise.',[samplePitch sigma sigma 100*noiseLevel]);
        logMessage('The magnitude of its Fourier transform abs(fft2(img)) ~ (%0.3f/f)exp(-f^2/(2*%0.3f^2)) + %0.6f.',[amplitude sigmaFt noiseLevelFt]);
    end
    
    
    
   
%     %%%%%%%%%%%%%%%%%%%
%     imgFtAbs=abs(fftshift(fft2(img)));
%     
%     nbAngles=7;
%     tfRange=[0:(nbAngles-1)]*pi/nbAngles;
%     nbRadii=sum(imgSize(1:2))/2;
%     rfRange=opticalCutOffFreq*([1:nbRadii]-1)./nbRadii;
%     [TFI,RFI]=ndgrid(tfRange,rfRange);
%     imgFtAbsPol=interp2(YF,XF,imgFtAbs,RFI.*cos(TFI),RFI.*sin(TFI),'*linear',0);
%     
%     
%     imgFtAbsPolNorm=imgFtAbsPol./noiseLevelFt;
%     imgFtAbsPolNoNoise=sqrt(max(0,imgFtAbsPolNorm.^2-1));
%         
%     MTF=RFI.*imgFtAbsPolNoNoise;
%     MTF=MTF./mean(MTF(:,2+10));
%     
% %     MTFFit=0*MTF;
% %     gFit=@(x) x(1)*exp(-abs(rfRange(2:end)).^2./(2*(1/(x(2)*2*pi)).^2));
% %     for tIdx=1:size(MTF,1),
% %         p=fminsearch(@(x) max(log(abs(gFit(x)-MTF(tIdx,2:end)))),[MTF(tIdx,2) 10],optimset('Display','iter'));
% %         p(2)
% % %         MTFFit(tIdx,:)=gFit(p);
% %     end
%     
%     
%     testSpectrum=fftshift(gaussian2DFt(xRange,yRange,4.5));
%     testSpectrum=testSpectrum./testSpectrum(1+floor(end/2),floor(end/2)+2+10);
%     testSpectrum=max(0.01,testSpectrum);
%     
%     close all;
%     figure;
%     semilogy(RFI.',MTF.'); hold on;
%     semilogy(RFI(1,:).',mean(MTF).','LineWidth',3,'Color',0.0*[1 1 1]); hold on;
%     semilogy(RF(1+floor(end/2),1+floor(end/2):end),testSpectrum(1+floor(end/2),1+floor(end/2):end),'LineWidth',3,'Color',[1 0 0]); hold on;
%     ylim([1e-2,10]);
%     
%     accurateFreqs=abs(XF)>4./imgSize(1)./samplePitch(1) & abs(YF)>4./imgSize(2)./samplePitch(2);
end

% function [GFt xfRange]=gaussianFt(xRange,sigma)
%     xfRange=calcFrequencyRanges(xRange);
%     sigmaFt=1./(sigma*2*pi);
%     GFt=exp(-abs(xfRange).^2./(2*sigmaFt^2));
% end

function [GFt xfRange yfRange]=gaussian2DFt(xRange,yRange,sigma)
    [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange);
    [XF,YF]=ndgrid(xfRange,yfRange);
    RF2=XF.^2+YF.^2;
    clear XF YF;
    sigmaFt=1./(sigma*2*pi);
    GFt=exp(-RF2./(2*sigmaFt^2));
end
