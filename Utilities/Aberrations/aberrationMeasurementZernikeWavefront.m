% [measuredPupilFunction zernikeCoefficients Rho Phi]=aberrationMeasurementZernikeWavefront(slm,selectedZernikeCoefficientIndexes,probeFunctor,progressFunctor)
%
% Determines the aberration by tweaking the low-order Zernike terms, hence it cannot correct for amplitude modulations
%
% Input:
%   slm: An SLM object
%   selectedZernikeCoefficientIndexes: the indexes of the standard Zernike polynomials to probe.
%   probeFunctor: a function that returns the probe value, it optionally takes as argument the complex field at each pixel of the SLM
%   progressFunctor: When specified, this function will be executed every spatial frequency probe, it takes as optional arguments:
%                    - fractionDone: the completion ratio
%                    - currentPupilFunctionEstimate: the complex matrix with the current estimate of the pupil function
%
% Returns:
%   measuredPupilFunction: the complex matrix with the estimate of the aberrated pupil function
%   zernikeCoefficients: the standard Zernik coefficients of the measured aberration.
%   Rho: the radial coordinate
%   Phi: the azimutal coordinate
%
% Example:
%     slmSize=[20 20]; % [rows cols]
%     referenceDeflectionFrequency=[1/4 1/4];
%     slm=PhaseSLM(2,referenceDeflectionFrequency,[0 0 slmSize]);
%     modulationFunctor=@(fieldModulationAtPupil) pause(1/30);
%     %Using:
%     function psf=calcPsf(fieldModulationAtPupil)
%         imgSize=size(fieldModulationAtPupil);
%         [X,Y]=meshgrid([1:imgSize(2)],[1:imgSize(1)]);
%         R2=((X-imgSize(2)/2).^2+(Y-imgSize(1)/2).^2)./(min(imgSize)/2/2).^2;
%         aperture=R2<1;
%         aperture=aperture.*exp(2*2i*pi*R2);
%         fieldPsf=fftshift(fft2(ifftshift(fieldModulationAtPupil.*aperture)));
%         psf=abs(fieldPsf).^2;
%     end
%     function value=simulateDetection(img,pos)
%         value=img(pos(1),pos(2))*(1+0.05*randn());
%     end
%     probeFunctor=@(fieldModulationAtPupil) simulateDetection(calcPsf(fieldModulationAtPupil),[5 5]);
%
%     selectedZernikeCoefficientIndexes=[1  2 3  4  5 6  7 8  9 10 11]; % piston,
%                  tip(x),tilt(y), defocus, astigmatism-diag,astigmatism-x,
%                  coma-y,coma-x,  trefoil-y,trefoil-x,  spherical
%
%     measuredPupilFunction=aberrationMeasurementZernikeWavefront(slm,selectedZernikeCoefficientIndexes,probeFunctor)
%
function [measuredPupilFunction zernikeCoefficients X Y Rho Phi]=aberrationMeasurementZernikeWavefront(slm,selectedZernikeCoefficientIndexes,probeFunctor,progressFunctor,referenceAberration)
    if (nargin<1)
        close all;
        figure();
        image(zeros(600,800));
        displayNumber=gca();
        
        referenceDeflectionFrequency=[1 1]/10;
        slm=PhaseSLM(displayNumber,referenceDeflectionFrequency);
        slm.stabilizationTime=0.01;
    end
    if (nargin<2)
        selectedZernikeCoefficientIndexes=[2 3 4]; % tip,tilt, and defocus only
    end
    slmSize=slm.regionOfInterest(3:4);
    if (nargin<3)
        probePos=1+floor(slmSize./2);
        probeFunctor=@(fieldModulationAtPupil) simulateDetection(calcPsf(fieldModulationAtPupil),probePos);
    end
    if (nargin<4)
        progressFunctor=@(fractionDone,currentPupilFunctionEstimate) progressFunction(fractionDone,currentPupilFunctionEstimate);
    end
    
    [X,Y]=meshgrid([1:slmSize(2)]-floor(slmSize(2)/2)-1,[1:slmSize(1)]-floor(slmSize(1)/2)-1);
    pupilRadius=sqrt(sum(ceil(slmSize./2).^2));
       
    if (nargin<5 || isempty(referenceAberration))
        referenceAberration=1;
    end
    if (isa(referenceAberration,'function_handle'))
        referenceAberration=referenceAberration(X-floor(slmSize(2)/2)-1,Y-floor(slmSize(1)/2)-1);
    end
    
    maxIterations=100;
    coefficientWeights=100;
    
    slm.correctionFunction=slm.correctionFunction.*referenceAberration;
    
    nbCoefficients=size(selectedZernikeCoefficientIndexes,2);
    selectedZernikeCoefficients0=zeros(1,nbCoefficients);
    
    function value=modulateAndProbe(selectedZernikeCoefficients)
        zernikeCoefficients=[];
        zernikeCoefficients(selectedZernikeCoefficientIndexes)=selectedZernikeCoefficients.*coefficientWeights;
        wavefront=zernikeComposition(X/pupilRadius,Y/pupilRadius,zernikeCoefficients);
        slm.modulate(exp(-2i*pi*wavefront));
        
        value = -probeFunctor(); % Maximize the intensity
    end

    function stop = outputFunction(selectedZernikeCoefficients,optimValues,state)
        zernikeCoefficients=[];
        zernikeCoefficients(selectedZernikeCoefficientIndexes)=selectedZernikeCoefficients.*coefficientWeights;
        zernikeCoefficients
        optimValues.fval
        %Report progress
        try
            switch (nargin(progressFunctor))
                case 0
                    cont=progressFunctor();
                case 1
                    cont=progressFunctor(optimValues.iteration/maxIterations);
                otherwise
                    wavefront=zernikeComposition(X/pupilRadius,Y/pupilRadius,zernikeCoefficients);
                    currentPupilFunctionEstimate=exp(2i*pi*wavefront);
                    cont=progressFunctor(optimValues.iteration/maxIterations,currentPupilFunctionEstimate);
            end
        catch TooManyLHSExc
            if (strcmp(TooManyLHSExc.identifier,'MATLAB:maxlhs') || strcmp(TooManyLHSExc.identifier,'MATLAB:TooManyOutputs'))
                cont=true; % continue until all probes are done or the progress functor returns false
            else
                % Error occured in the progress display, exiting...
                rethrow(TooManyLHSExc);
            end
        end
        
        stop=(cont==false); % Add some fool-proving for incorrect progressFunctor's
    end

    % Optimize now
    selectedZernikeCoefficients=fminsearch(@modulateAndProbe,selectedZernikeCoefficients0,optimset('MaxIter',maxIterations,'TolX',1e-6,'Display','iter','OutputFcn',@outputFunction));
    zernikeCoefficients=[];
    zernikeCoefficients(selectedZernikeCoefficientIndexes)=selectedZernikeCoefficients.*coefficientWeights;
    wavefront=zernikeComposition(X/pupilRadius,Y/pupilRadius,zernikeCoefficients);
    measuredPupilFunction=exp(2i*pi*wavefront);
    
    if (nargout==0)
        showImage(measuredPupilFunction,-1);
        clear('measuredPupilFunction');
    end
    if (nargout>4)
        [Phi Rho]=cart2pol(X/pupilRadius,Y/pupilRadius);
    end
end

function progressFunction(fractionDone,currentPupilFunctionEstimate)
    global debugOutput;
    
    persistent prevPctDone;
    pctDone=floor(100*fractionDone);
    if (~exist('prevPctDone') || isempty(prevPctDone) || prevPctDone~=pctDone)
        logMessage('%u%% done.',pctDone);
        prevPctDone=pctDone;
    end
    
    if (debugOutput)
        currentPupilFunctionEstimate=currentPupilFunctionEstimate./max(abs(currentPupilFunctionEstimate(:)));
        SNR=100;
        pupilFunctionCorrection=conj(currentPupilFunctionEstimate)./(abs(currentPupilFunctionEstimate).^2+(1/SNR).^2);
        nbGrayLevels=20;
        pupilFunctionCorrection=exp(1i*angle(pupilFunctionCorrection)).*floor(nbGrayLevels*min(1,abs(pupilFunctionCorrection)))./nbGrayLevels;
        psf=calcPsf(pupilFunctionCorrection);
        subplot(2,2,2);
        showImage(psf./max(abs(psf(:))),[],[],[],gca);
        drawnow();
    end
end

function [value auxValues]=displayOnSLMAndAcquireSignal(slm,sampleProbe,probeFunctor,referenceProbe)
    %Display mask
    combinedProbes=referenceProbe+sampleProbe;
    slm.modulate(combinedProbes);
    
    switch (nargout(probeFunctor))
        case 1
            if (nargin(probeFunctor)==0)
                value=probeFunctor();
            else
                value=probeFunctor(combinedProbes);
            end
        otherwise
            if (nargin(probeFunctor)==0)
                [value auxValues]=probeFunctor();
            else
                [value auxValues]=probeFunctor(combinedProbes);
            end
    end
end

function psf=calcPsf(fieldModulationAtPupil)
    %Simulate defocus
    w20=3;
    imgSize=size(fieldModulationAtPupil);
    
    [X,Y]=meshgrid([1:imgSize(2)],[1:imgSize(1)]);
    R2=((X-imgSize(2)/2).^2+(Y-imgSize(1)/2).^2)./(min(imgSize)/2).^2;
    aperture=R2<(1/2)^2;
    openFraction=sum(aperture(:))/numel(aperture);
    aperture=aperture.*exp(2i*pi*w20*R2);
    
    fieldPsf=fftshift(fft2(ifftshift(fieldModulationAtPupil.*aperture)))/sqrt(numel(fieldModulationAtPupil)^2*openFraction);
    psf=abs(fieldPsf).^2; % Integral 1
    
    psf=psf*1000;
end

function [value auxValues]=simulateDetection(img,pos)
    img=img.*(1+0.05*randn(size(img))); %Simulate photon noise
    value=img(pos(1),pos(2));
    auxValues=img;
end
