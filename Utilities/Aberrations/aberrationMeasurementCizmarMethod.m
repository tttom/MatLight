% [measuredPupilFunction pupilFunctionEstimates probeField4DMatrix samplePosX samplePosY]=aberrationMeasurementCizmarMethod(slm,probeGridSize,probeFunctor,progressFunctor)
%
% Determines the aberration using Tomas Cizmar's method as published in Nature Photonics.
%
% Input:
%   slm: An Slm object
%   probeGridSize: a vector of three scalars, the number of probes in y and x followed by the number of phase probes (must be >=3)
%   probeFunctor: a function that returns the probe value, it optionally takes as argument the complex field at each pixel of the SLM
%   progressFunctor: When specified, this function will be executed every spatial frequency probe, it takes as optional arguments:
%                    - fractionDone: the completion ratio
%                    - currentPupilFunctionEstimate: the complex matrix with the current estimate of the pupil function
%
% Returns:
%   measuredPupilFunction: the complex matrix with the estimate of the aberrated pupil function
%   pupilFunctionEstimates: the eigenvector maximizing the probe value, ordered into a matrix where the rows and columns correspond to the spatial frequencies of the probes as specified by samplePosX and samplePosY.
%   probeField4DMatrix: the same as pupilFunctionEstimates, but for every auxilary value returned by the probeFunctor. If the auxilary value is a matrix, it will be stacked along the higher dimensions of this matrix. Typical usage would be that the probeFunctor returns the full image.
%   samplePosX: the horizontal spatial frequency in cycles/SLM pixel for each probe
%   samplePosY: the vertical spatial frequency in cycles/SLM pixel for each probe
%   referenceProbeAmplitude: The amplitude of the reference beam at the probe position. This is required if the reference caries an additional aberration.
%   referenceAuxProbesAmplitude: The amplitude of the reference beam at the auxilary reference probe positions. This is required if the reference caries an additional aberration.
%
%
% Example:
%     slmSize=[20 20]; % [rows cols]
%     referenceDeflectionFrequency=[1/4 1/4];
%     slm=PhaseSLM(2,referenceDeflectionFrequency,[0 0 slmSize]);
%     nbOfPhaseProbes=3;
%     probeGridSize=[20 20 nbOfPhaseProbes];
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
%     probeFunctor=@(fieldModulationAtPupil) simulateDetection(calcPsf(fieldModulationAtPupil),[15,15]);
%
%     measuredPupilFunction=aberrationMeasurementCizmarMethod(slm,probeGridSize,probeFunctor)
%
function [measuredPupilFunction pupilFunctionEstimates probeField4DMatrix samplePosX samplePosY referenceProbeAmplitude referenceAuxProbesAmplitude]=aberrationMeasurementCizmarMethod(slm,probeGridSize,probeFunctor,progressFunctor,referenceAberration)
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
        probeGridSize=[25 25 3]; %[12 16 3]; % [rows,cols,phases]
    end
    if (length(probeGridSize)<2)
        probeGridSize(2)=0;
        logMessage('The horizontal deflection has not been specified, defaulting to 0 (vertical only deflection).');
    end
    if (length(probeGridSize)<3)
        probeGridSize(3)=3;
        logMessage('The number of phases to sample has not been specified, defaulting to %u',probeGridSize(3));
    end
    slmSize=slm.regionOfInterest(3:4);
    if (nargin<3)
        probePos=1+floor(slmSize./2);
        probeFunctor=@(fieldModulationAtPupil) simulateDetection(calcPsf(fieldModulationAtPupil),probePos);
    end
    if (nargin<4)
        progressFunctor=@(fractionDone,currentPupilFunctionEstimate) progressFunction(fractionDone,currentPupilFunctionEstimate);
    end
    
    nbOfPhaseProbes=probeGridSize(3);
    probeSpacing=floor(slmSize./probeGridSize(1:2)); % [y, x] = [row,col]
    probeSize=probeSpacing;
    probeShapeFunctor=@(X,Y) double(X>=-probeSize(1)/2 & X<probeSize(1)/2 & Y>=-probeSize(2)/2 & Y<probeSize(2)/2);
    [X,Y]=meshgrid([1:slmSize(2)]-floor(slmSize(2)/2)-1,[1:slmSize(1)]-floor(slmSize(1)/2)-1);
    
    [samplePosX,samplePosY]=meshgrid(probeSpacing(2)*([0.5:probeGridSize(2)]-probeGridSize(2)/2),probeSpacing(1)*([0.5:probeGridSize(1)]-probeGridSize(1)/2));
    sampleR2=samplePosX.^2+samplePosY.^2;
    [ign sortI]=sort(sampleR2(:));
        
    if (nargin<5 || isempty(referenceAberration))
        referenceAberration=1;
    end
    if (isa(referenceAberration,'function_handle'))
        referenceAberration=referenceAberration(X-floor(slmSize(2)/2)-1,Y-floor(slmSize(1)/2)-1);
    end
    
    slm.correctionFunction=slm.correctionFunction.*referenceAberration;
    
    referenceProbe=probeShapeFunctor(X-samplePosX(sortI(1)),Y-samplePosY(sortI(1)));
    
    [ign randomizedPhaseI]=sort(rand(1,nbOfPhaseProbes));
    
    %Probe for different deflections
    currentPupilFunctionEstimate=ones(slmSize);
    pupilFunctionEstimates=ones(size(samplePosX));
    probeField4DMatrix=[];
    auxProbesCoefficientInMatrixForm=[];
    nbSampleDeflections=prod(probeGridSize(1:2));
    referenceProbeAmplitude=[];
    referenceAuxProbesAmplitude=[];
    for sampleProbeIdx=1:nbSampleDeflections,
        sampleProbe=probeShapeFunctor(X-samplePosX(sortI(sampleProbeIdx)),Y-samplePosY(sortI(sampleProbeIdx)));
        coincidenceFactor=1+double(any(abs(referenceProbe(:)+sampleProbe(:))>1));
        %Test different phase offsets
        newValues=zeros(1,nbOfPhaseProbes);
        auxValues={};
        for phaseIdx=randomizedPhaseI,
            if (nargout(probeFunctor)==1 || nargout<3)
                newValues(phaseIdx)=displayOnSLMAndAcquireSignal(slm,sampleProbe*exp(2i*pi*(phaseIdx-1)/nbOfPhaseProbes)/coincidenceFactor,probeFunctor,referenceProbe/coincidenceFactor);
            else
                [newValues(phaseIdx) auxValues{phaseIdx}]=displayOnSLMAndAcquireSignal(slm,sampleProbe*exp(2i*pi*(phaseIdx-1)/nbOfPhaseProbes)/coincidenceFactor,probeFunctor,referenceProbe/coincidenceFactor);
                auxValues{phaseIdx}=(coincidenceFactor^2)*auxValues{phaseIdx};
            end
        end
        newValues=(coincidenceFactor^2)*newValues;
        
        %Work out the phase and amplitude from the sampling using 'lock-in' amplification
        probeCoefficient=newValues*exp(-2i*pi*([1:nbOfPhaseProbes]-1)/nbOfPhaseProbes).';
        %Do the same for the auxilary probes
        if (nargout>2 && ~isempty(auxValues))
            auxProbesCoefficientInMatrixForm=auxValues{1};
            for phaseIdx=2:nbOfPhaseProbes,
                % Correct for reduction in referenceFraction
                auxProbesCoefficientInMatrixForm=auxProbesCoefficientInMatrixForm+auxValues{phaseIdx}*exp(-2i*pi*(phaseIdx-1)/nbOfPhaseProbes);
            end
        end
        
        %Verify and store the reference, this is always calculated first
        if (isempty(referenceProbeAmplitude))
            %The first probe is a self-reference test, so the result should be a positive real number
            probeErrorEstimate=abs(angle(probeCoefficient));
            if (probeErrorEstimate>0.01)
                logMessage('The estimated measurement error is large: %0.2f%%',probeErrorEstimate*100);
            end
            if (~isempty(auxValues))
                auxMatrixErrorEstimate=sqrt(mean(abs(angle(auxProbesCoefficientInMatrixForm(:))).^2));
                if (auxMatrixErrorEstimate>0.01)
                    logMessage('The estimated measurement error for the auxilary probes is large: %0.2f%%',auxMatrixErrorEstimate*100);
                end
            end
            
            %Force to be real, the imaginary must be a measurement error
            probeCoefficient=real(probeCoefficient);
            if (~isempty(auxProbesCoefficientInMatrixForm))
                auxProbesCoefficientInMatrixForm=real(auxProbesCoefficientInMatrixForm);
            end
            
            %Store the reference
            referenceProbeAmplitude=sqrt(max(0,probeCoefficient));
            if (~isempty(auxValues))
                referenceAuxProbesAmplitude=sqrt(max(0,auxProbesCoefficientInMatrixForm));
            end
        end
        
        %Store the measurements
        if (sampleProbeIdx>0)
            %The probe scalar
            pupilFunctionEstimates(sortI(sampleProbeIdx))=probeCoefficient;
            currentPupilFunctionEstimate=currentPupilFunctionEstimate+probeCoefficient*conj(sampleProbe);
            %Store any auxilary values as well
            if (~isempty(auxProbesCoefficientInMatrixForm))
                if (isempty(probeField4DMatrix))
                    probeField4DMatrix=zeros([size(samplePosX) size(auxProbesCoefficientInMatrixForm)]);
                end
                probeField4DMatrix(sortI(sampleProbeIdx)+numel(samplePosX)*[0:numel(auxProbesCoefficientInMatrixForm)-1])=auxProbesCoefficientInMatrixForm;
            end
                
            %Report progress
            try
                switch (nargin(progressFunctor))
                    case 0
                        cont=progressFunctor();
                    case 1
                        cont=progressFunctor(sampleProbeIdx/nbSampleDeflections);
                    otherwise
                        cont=progressFunctor(sampleProbeIdx/nbSampleDeflections,currentPupilFunctionEstimate);
                end
            catch TooManyLHSExc
                if (strcmp(TooManyLHSExc.identifier,'MATLAB:maxlhs') || strcmp(TooManyLHSExc.identifier,'MATLAB:TooManyOutputs'))
                    cont=true; % continue until all probes are done or the progress functor returns false
                else
                    % Error occured in the progress display, exiting...
                    rethrow(TooManyLHSExc);
                end
            end
        end
        cont=(cont~=false); % Add some fool-proving for incorrect progressFunctor's
        if (~cont)
            break;
        end
    end % Measurement done
    
    %Normalize the measurements to a maximum of unity transmission
    pupilFunctionEstimates(sortI(1:sampleProbeIdx))=pupilFunctionEstimates(sortI(1:sampleProbeIdx))./max(abs(pupilFunctionEstimates(sortI(1:sampleProbeIdx))));
    
    % Prepare output
%     measuredPupilFunction=currentPupilFunctionEstimate; %Will not work if the probe size is not equal to the probe spacing
    % Interpolate amplitude and phase separately
    measuredPupilFunctionAbs=interp2(samplePosX,samplePosY,abs(pupilFunctionEstimates),X,Y,'*nearest',0);
    measuredPupilFunction=interp2(samplePosX,samplePosY,pupilFunctionEstimates,X,Y,'*nearest',0);
    measuredPupilFunction=measuredPupilFunctionAbs.*exp(1i*angle(measuredPupilFunction));
    
    %If no special reference aberration has been given, it is that of the zero deflection
    if (isempty(referenceProbeAmplitude))
        referenceProbe=max(0,real(pupilFunctionEstimates(floor(end/2)+1,floor(end/2)+1)));
        referenceProbeAmplitude=sqrt(referenceProbe);
        if (~isempty(probeField4DMatrix))
            referenceAuxProbes=max(0,real(probeField4DMatrix(floor(end/2)+1,floor(end/2)+1,:,:))); %size == [1 1 roiSize(1:2)]
            referenceAuxProbesAmplitude=sqrt(squeeze(referenceAuxProbes));
        end
    end
    
    if (nargout==0)
        showImage(measuredPupilFunction,-1);
        clear('measuredPupilFunction');
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
