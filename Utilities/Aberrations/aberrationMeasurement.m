% [measuredPupilFunction eigenVector probeField4DMatrix sampleX sampleY]=aberrationMeasurement(slm,probeGridSize,probeFunctor,progressFunctor)
%
% Determines the aberration based on interference of different deflections
% and calculates the correction function. Multiple probes can be used, such
% as individual pixels of a camera, or NSOM tip. See
% camAberrationMeasurement.m for an example. When multiple probes are
% given, the aberrations are estimated simulatanously for each probe.
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
%   eigenVector: the eigenvector maximizing the probe value, ordered into a matrix where the rows and columns correspond to the spatial frequencies of the probes as specified by sampleX and sampleY.
%   probeField4DMatrix: the same as eigenVector, but for every auxilary value returned by the probeFunctor. If the auxilary value is a matrix, it will be stacked along the higher dimensions of this matrix. Typical usage would be that the probeFunctor returns the full image.
%   sampleX: the horizontal spatial frequency in cycles/SLM pixel for each probe
%   sampleY: the vertical spatial frequency in cycles/SLM pixel for each probe
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
%     measuredPupilFunction=aberrationMeasurement(slm,probeGridSize,probeFunctor)
%
function [measuredPupilFunction eigenVector probeField4DMatrix sampleX sampleY referenceProbeAmplitude referenceAuxProbesAmplitude]=aberrationMeasurement(slm,probeGridSize,probeFunctor,progressFunctor,referenceAberration)
    if (nargin<1)
        referenceDeflectionFrequency=[1 -1]*1/10;
        slm=PhaseSLM(2,referenceDeflectionFrequency,[0 0 60 40]);
    end
    if (nargin<2)
        probeGridSize=[15 15 3]; % [rows,cols,phases]
    end
    if (length(probeGridSize)<2)
        probeGridSize(2)=0;
        logMessage('The horizontal deflection has not been specified, defaulting to 0 (vertical only deflection).');
    end
    if (length(probeGridSize)<3)
        probeGridSize(3)=3;
        logMessage('The number of phases to sample has not been specified, defaulting to %u',probeGridSize(3));
    end
    slmSize=slm.regionOfInterest;
    slmSize=slmSize(3:4);
    if (nargin<3)
        probePos=1+floor(slmSize./2)+round(slm.referenceDeflectionFrequency(1:2).*slmSize)+[0 -5];
        probeFunctor=@(fieldModulationAtPupil) simulateDetection(calcPsf(fieldModulationAtPupil),probePos);
    end
    
    %the base deflection for a phase SLM in cycles per pixel, [cycles/row cycles/col]
    samplingDeflectionFrequencyStep=[1 1]./slmSize; % [y, x] = [row,col]
    
    nbOfPhaseProbes=probeGridSize(3);

    [X,Y]=meshgrid([1:slmSize(2)],[1:slmSize(1)]);
    if (nargin<4)
        progressFunctor=@(fractionDone,currentPupilFunctionEstimate) progressFunction(fractionDone,currentPupilFunctionEstimate);
    end
    [sampleX,sampleY]=meshgrid([-floor(probeGridSize(2)/2):floor((probeGridSize(2)-1)/2)]*samplingDeflectionFrequencyStep(2),[-floor(probeGridSize(1)/2):floor((probeGridSize(1)-1)/2)]*samplingDeflectionFrequencyStep(1));
    sampleR2=sampleX.^2+sampleY.^2;
    [ign sortI]=sort(sampleR2(:));
        
    if (nargin<5 || isempty(referenceAberration))
        referenceAberration=1;
    end
    if (isa(referenceAberration,'function_handle'))
        referenceAberration=referenceAberration(X-floor(slmSize(2)/2)-1,Y-floor(slmSize(1)/2)-1);
    end

%     %Pick exactly 50% of the pixels
%     slmPixelsForReferenceField=mod(X+Y,2); %checker-board pattern encoding
    %Start with half the field in the reference beam
    referenceFraction=0.5;
    
    % Randomize phase probes to limit a bias due to laser fluctuations
    [ign randomizedPhaseI]=sort(rand(1,nbOfPhaseProbes));
    
    %Probe for different deflections
    currentPupilFunctionEstimate=zeros(slmSize);
    eigenVector=zeros(size(sampleX));
    probeField4DMatrix=[];
    auxProbesCoefficientInMatrixForm=[];
    nbSampleDeflections=prod(probeGridSize(1:2));
    referenceProbeAmplitude=[];
    referenceAuxProbesAmplitude=[];
    for sampleDeflectionIdx=double(all(referenceAberration(:)==1)):nbSampleDeflections,
        if (sampleDeflectionIdx>0)
            samplingDeflectionFrequency=[sampleY(sortI(sampleDeflectionIdx)) sampleX(sortI(sampleDeflectionIdx))];
            differentialDeflection=exp(2i*pi*(X*samplingDeflectionFrequency(2)+Y*samplingDeflectionFrequency(1)));
        else
            %If an aberration of the reference is specified, probe this first
            differentialDeflection=referenceAberration;
        end
        %Test different phase offsets
        newValues=zeros(1,nbOfPhaseProbes);
        auxValues={};
        for phaseIdx=randomizedPhaseI,
            if (nargout(probeFunctor)==1 || nargout<3)
                newValues(phaseIdx)=displayOnSLMAndAcquireSignal(referenceFraction,slm,differentialDeflection*exp(2i*pi*(phaseIdx-1)/nbOfPhaseProbes),probeFunctor,referenceAberration);
            else
                [newValues(phaseIdx) auxValues{phaseIdx}]=displayOnSLMAndAcquireSignal(referenceFraction,slm,differentialDeflection*exp(2i*pi*(phaseIdx-1)/nbOfPhaseProbes),probeFunctor,referenceAberration);
            end
        end
        
        referenceFractionCorrection=1/(referenceFraction*(1-referenceFraction));
        
        %Work out the phase and amplitude from the sampling using 'lock-in' amplification
        probeCoefficient=referenceFractionCorrection*newValues*exp(-2i*pi*([1:nbOfPhaseProbes]-1)/nbOfPhaseProbes).';
%         probeScaling=mean(newValues); % proportional with rrAA+(1-r)(1-r)BB
        % maximize  r*(1-r)/sqrt(r^2+((1-r)*x)^2) for r=refFrac, with x is the fraction B/A
        
        if (nargout>2 && ~isempty(auxValues))
            auxProbesCoefficientInMatrixForm=auxValues{1};
            for phaseIdx=2:nbOfPhaseProbes,
                % Correct for reduction in referenceFraction
                auxProbesCoefficientInMatrixForm=auxProbesCoefficientInMatrixForm+referenceFractionCorrection*auxValues{phaseIdx}*exp(-2i*pi*(phaseIdx-1)/nbOfPhaseProbes);
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
        if (sampleDeflectionIdx>0)
            %The probe scalar
            eigenVector(sortI(sampleDeflectionIdx))=probeCoefficient;
            currentPupilFunctionEstimate=currentPupilFunctionEstimate+probeCoefficient*conj(differentialDeflection);
            %Store any auxilary values as well
            if (~isempty(auxProbesCoefficientInMatrixForm))
                if (isempty(probeField4DMatrix))
                    probeField4DMatrix=zeros([size(sampleX) size(auxProbesCoefficientInMatrixForm)]);
                end
                probeField4DMatrix(sortI(sampleDeflectionIdx)+numel(sampleX)*[0:numel(auxProbesCoefficientInMatrixForm)-1])=auxProbesCoefficientInMatrixForm;
            end
                
            %Report progress
            try
                switch (nargin(progressFunctor))
                    case 0
                        cont=progressFunctor();
                    case 1
                        cont=progressFunctor(sampleDeflectionIdx/nbSampleDeflections);
                    otherwise
                        cont=progressFunctor(sampleDeflectionIdx/nbSampleDeflections,currentPupilFunctionEstimate);
                end
            catch TooManyLHSExc
                if (strcmp(TooManyLHSExc.identifier,'MATLAB:maxlhs'))
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
    end
    
    measuredPupilFunction=currentPupilFunctionEstimate; % or calcPupilFunctionFromProbes(slmSize,eigenVector);
    
%     %Re-interpolate the pixels on the SLM that were used for the reference beam (absolute value and argument separately to avoid a bias towards lower absolute amplitudes)
%     %for measuredPupilFunction:
%     measuredPupilFunctionAbs=abs(measuredPupilFunction);
%     interpolatedFunctionAbs=(measuredPupilFunctionAbs([1 1:end-1],:)+measuredPupilFunctionAbs([2:end end],:)+measuredPupilFunctionAbs(:,[1 1:end-1])+measuredPupilFunctionAbs(:,[2:end end]))./4;
%     interpolatedFunctionArg=angle(measuredPupilFunction([1 1:end-1],:)+measuredPupilFunction([2:end end],:)+measuredPupilFunction(:,[1 1:end-1])+measuredPupilFunction(:,[2:end end]));
%     interpolatedFunction=exp(1i*interpolatedFunctionArg).*interpolatedFunctionAbs;
%     measuredPupilFunction(slmPixelsForReferenceField>0)=interpolatedFunction(slmPixelsForReferenceField>0);
    
    %If no special reference aberration has been given, it is that of the zero deflection
    if (isempty(referenceProbeAmplitude))
        referenceProbe=max(0,real(eigenVector(floor(end/2)+1,floor(end/2)+1)));
        referenceProbeAmplitude=sqrt(referenceProbe);
        if (~isempty(probeField4DMatrix))
            referenceAuxProbes=max(0,real(probeField4DMatrix(floor(end/2)+1,floor(end/2)+1,:,:))); %size == [1 1 roiSize(1:2)]
            referenceAuxProbesAmplitude=sqrt(squeeze(referenceAuxProbes));
        end
    end
    
    if (nargout==0)
        showImage(measuredPupilFunction./max(abs(measuredPupilFunction(:))));
        clear('measuredPupilFunction');
    end
end

function cont=progressFunction(fractionDone,currentPupilFunctionEstimate)
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
    end
    
    cont=true;
end

function [value auxValues]=displayOnSLMAndAcquireSignal(referenceFraction,slm,sampleDeflection,probeFunctor,referenceAberration)
    %Display mask
%     combinedDeflection=referenceAberration.*slmPixelsForReferenceField + sampleDeflection.*(1-slmPixelsForReferenceField); %checker-board pattern encoding
    combinedDeflection=referenceFraction*referenceAberration + (1-referenceFraction)*sampleDeflection;
    slm.modulate(combinedDeflection);
    
    switch (nargout(probeFunctor))
        case 1
            if (nargin(probeFunctor)==0)
                value=probeFunctor();
            else
                value=probeFunctor(combinedDeflection);
            end
        otherwise
            if (nargin(probeFunctor)==0)
                [value auxValues]=probeFunctor();
            else
                [value auxValues]=probeFunctor(combinedDeflection);
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
    aperture=aperture.*exp(2i*pi*w20*R2);
%     aperture=aperture.*exp(-R2./(.25.^2));
%     aperture=aperture.*R2; %sinc(sqrt(R2)./.25);
    fieldPsf=fftshift(fft2(ifftshift(fieldModulationAtPupil.*aperture)));
    psf=abs(fieldPsf).^2;
end

function [value auxValues]=simulateDetection(img,pos)
    img=img.*(1+0.05*randn(size(img))); %Simulate photon noise
    value=img(pos(1),pos(2));
    auxValues=img;
end
