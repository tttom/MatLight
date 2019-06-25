function gerschbergSaxton()
    slm={};
    slm.minGerschbergSaxtonIterations=10;
    slm.maxGerschbergSaxtonIterations=25;
    slm.maxGerschbergSaxtonError=0.01;
    slm.twoPiEquivalent=1.5;
    slm.crossTalkStds=[1 1]*1.5;
    slm.equalPhaseBlockSize=[1 1]*1;
    
    
    gridSize=[600 800];
    slm.referenceDeflectionFrequency=[1/10 1/10];
    xRange=([1:gridSize(1)]-1-floor(gridSize(1)/2));
    yRange=([1:gridSize(2)]-1-floor(gridSize(2)/2));
    [XPupil,YPupil]=ndgrid(xRange,yRange);
    slm.referenceDeflectionAngle=2*pi*(XPupil*slm.referenceDeflectionFrequency(2)+YPupil*slm.referenceDeflectionFrequency(1));
%     targetPupil=parsePupilEquation('0.3*(R>.9 & R<1)',2*XPupil/min(gridSize),2*YPupil/min(gridSize));
    targetPupil=parsePupilEquation('0.25*(R<1 & R>0.5)*exp(1i*P)',2*XPupil/min(gridSize),2*YPupil/min(gridSize));

    tic;
    initialModulation=checkerModulate(targetPupil,slm.referenceDeflectionAngle,slm.equalPhaseBlockSize);
    [pupilModulation err nbIterations targetField optimizedField]=optimizePhaseOnlyPupilFunction(slm,initialModulation,targetPupil);
    toc
    pupil=exp(2i*pi*pupilModulation/slm.twoPiEquivalent);
    
    normalizationForDisplay=.90./max(abs(targetField(:)));

    close all;
    figure;
    axs(1)=subplot(2,2,1);
    showImage(targetPupil+1e-6i);
    axs(3)=subplot(2,2,3);
    showImage(targetField*normalizationForDisplay+1e-6i);
    axs(2)=subplot(2,2,2);
    showImage(pupil);
    axs(4)=subplot(2,2,4);
    showImage(optimizedField*normalizationForDisplay+1e-6i);
    linkaxes(axs([1 2]));
    linkaxes(axs([3 4]));
    drawnow();
    
end

% Uses the Gerschberg-Saxton algorithm to optimize the pupil function to
% minimize the L2-difference between the image of the phase only pupil
% function and that of the target pupil function.
function [pupilFunctionPhase err nbIterations targetField optimizedField]=optimizePhaseOnlyPupilFunction(slm,initialPupilModulation,targetPupilFunction)
    pupilModulation=initialPupilModulation;
    clear initialPupilModulation;
    targetField=propagateForward(targetPupilFunction,slm.referenceDeflectionFrequency);
    clear targetPupil;
    
    % Iterate towards solution
    nbIterations=0; err=Inf;
    pupilFunction=exp(2i*pi*pupilModulation);
    while (nbIterations<min(slm.minGerschbergSaxtonIterations,slm.maxGerschbergSaxtonIterations) || (err>slm.maxGerschbergSaxtonError && nbIterations<slm.maxGerschbergSaxtonIterations))
        nbIterations=nbIterations+1;
        %Restrict, modulate, and propagate
        restrictedPupilField=calcPupilFieldFromModulation(pupilModulation,slm.referenceDeflectionAngle,slm.twoPiEquivalent,slm.crossTalkStds);
        optimizedField=propagateForward(restrictedPupilField,slm.referenceDeflectionFrequency);
        %Correct
        differenceInImage=optimizedField-targetField;
        differenceInPupil=propagateBackward(differenceInImage,size(pupilModulation));
        pupilFunction=pupilFunction-differenceInPupil;
        pupilModulation=angle(pupilFunction)/(2*pi);
        
        err=norm(differenceInImage(:));
        [nbIterations err]
    end
    pupilFunctionPhase=mod(angle(pupilFunction)/(2*pi),1);
end

function outField=propagateForward(inField,imageSizeFraction)
    imageFractionSize=ceil(size(inField).*imageSizeFraction);
    xRange=[0:imageFractionSize(1)-1]-floor(imageFractionSize(1)/2);
    yRange=[0:imageFractionSize(2)-1]-floor(imageFractionSize(2)/2);
    inputSize=size(inField);
    inField(2*end,2*end)=0; inField=circshift(inField,inputSize-floor(inputSize/2)); %Zero pad and recenter
    outField=fftshift(ifft2(ifftshift(inField)));
    outField=outField(1+floor(end/2)+xRange,1+floor(end/2)+yRange);
    % or
%     outField=conj(czt2fromRanges(conj(inField),xRange,yRange))/sqrt(prod(size(inField)));
end
function outField=propagateBackward(inField,slmSize)
    inputSize=size(inField);
    if (any(inputSize<slmSize*2))
        inField(slmSize(1)*2,slmSize(2)*2)=0;
        inField=circshift(inField,slmSize-floor(inputSize/2));
    else
        if (any(inputSize>slmSize*2))
            inField=circshift(inField,slmSize-floor(inputSize/2));
            inField=inField(1:slmSize(1)*2,1:slmSize(2)*2);
        end
    end
    outField=circshift(fft2(ifftshift(inField)),floor(slmSize/2));
    outField=outField(1:end/2,1:end/2);
    % or
%     xRange=([0:slmSize(1)-1]-floor(slmSize(1)/2))*floor(inputSize(1)/2)/floor(slmSize(1)/2);
%     yRange=([0:slmSize(2)-1]-floor(slmSize(2)/2))*floor(inputSize(2)/2)/floor(slmSize(2)/2);
%     outField=czt2fromRanges(inField,xRange,yRange)*sqrt(prod(slmSize));
end

function restrictedInField=calcPupilFieldFromModulation(pupilModulation,referenceDeflectionAngle,twoPiEquivalent,crossTalkStds)
    pupilField=exp(2i*pi*pupilModulation);
    %Restrict phase
    maxPhase=pi/twoPiEquivalent;
    phase=mod(angle(pupilField)+referenceDeflectionAngle+pi,2*pi)-pi;
    phase(phase>maxPhase)=maxPhase;
    phase(phase<-maxPhase)=-maxPhase;
%     %Simulate cross talk
    crossTalkStds=abs(crossTalkStds);
    if (any(crossTalkStds>0))
        if (all(crossTalkStds>0))
            phase=convolve(phase,@(X,Y) exp(-0.5*((X/crossTalkStds(1)).^2+(Y/crossTalkStds(2)).^2)));
        else
            if (crossTalkStds(1)>0)
                phase=convolve(phase,@(X,Y) exp(-0.5*((X/crossTalkStds(1)).^2)));
            else
                phase=convolve(phase,@(X,Y) exp(-0.5*((Y/crossTalkStds(2)).^2)));
            end
        end
    end
    %Restrict amplitude
    restrictedInField=exp(1i*(phase-referenceDeflectionAngle));
end

function convolved=convolve(img,kernelFunctor)
    imgSize=size(img);
    [X,Y]=ndgrid([1:imgSize(1)]-floor(imgSize(1)/2)-1,[1:imgSize(2)]-floor(imgSize(2)/2)-1);
    
    imgFFT=fft2(img);
    kernelFFT=fft2(ifftshift(kernelFunctor(X,Y)));
    kernelFFT=kernelFFT./kernelFFT(1);
    convolved=ifft2(imgFFT.*kernelFFT,'symmetric');
end

function pupilModulation=checkerModulate(targetPupil,referenceDeflectionAngle,equalPhaseBlockSize)
    pupilSize=size(targetPupil);
    selectedPixels=mod(...
        (mod(repmat([1:pupilSize(1)].'-1,[1 pupilSize(2)])/equalPhaseBlockSize(1),2)>=1)+ ...
        (mod(repmat([1:pupilSize(2)]-1,[pupilSize(1) 1])/equalPhaseBlockSize(2),2)>=1),...
        2)>0;
    
    phaseValues=angle(targetPupil);
    phaseValues=mod((phaseValues+referenceDeflectionAngle)./(2*pi),1); %Convert to pixel value [0 1)
    amplitudeValues=abs(targetPupil);
    %make sure that the amplitude is in the dynamic range of the SLM
    if (max(abs(amplitudeValues(:)))>1)
        if (max(abs(amplitudeValues(:)))-1 > 4*eps(1))
            %If not a rounding error, signal this.
            logMessage('Warning: amplitude above 1, clipping it!');
        end
        amplitudeValues=min(amplitudeValues,1); %Convert to pixel value [0 1]
    end
    phaseDeviationFromMeanForAmplitude=acos(amplitudeValues)./(2*pi); %this is still the absolute value
    phaseDeviationFromMeanForAmplitude=phaseDeviationFromMeanForAmplitude.*(2*selectedPixels-1);
    imageForSLM=mod(phaseValues+phaseDeviationFromMeanForAmplitude,1);
    pupilModulation=mod(imageForSLM-referenceDeflectionAngle/(2*pi)+.5,1)-.5;
end
