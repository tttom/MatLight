% pupilFunctionCorrection=calcCorrectionFromPupilFunction(measuredPupilFunction,amplificationLimit)
%
% Calculates the complex value that the SLM has to modulate to counteract
% the aberration specified in measuredPupilFunction. The argument will thus
% be the inverse, and the amplitude will approximatelly be the reciprocal,
% however the signal reduction is limited by amplificationLimit to avoid
% blocking large parts of the pupil when it is illuminated unevenly.
%
% measuredPupilFunction: a complex matrix with the to-be-corrected
%     aberration at each SLM pixel
% amplificationLimit: The maximum dynamic range of the amplitude
%     modulation. An N-fold amplification of some parts of the pupil means
%     that part of the pupil has to be suppressed N-fold. To avoid low
%     signal this can be limited. (default: 10)
function pupilFunctionCorrection=calcCorrectionFromPupilFunction(measuredPupilFunction,amplificationLimit)
    if (nargin<2 || isempty(amplificationLimit))
        amplificationLimit=10;
    end

    measuredPupilFunction=amplificationLimit*measuredPupilFunction./max(abs(measuredPupilFunction(:)));
    
    pupilAmplitude=max(1,abs(measuredPupilFunction));
    pupilPhase=angle(measuredPupilFunction);
    
    pupilFunctionCorrection=exp(-1i*pupilPhase)./pupilAmplitude;
end