% [spectralWidthMeasuredInWavelength spectralWidthMeasuredInHz] = calcGaussianSpectrumWidthInWavelengthFromPulseDuration(pulseDurationMeasured, centerWavelength, propagationSpeed)
%
% Determines the bandwidth for Gaussian pulses of a certain duration.
%
% Inputs:
%     pulseDurationMeasured: the length of the pulse in time, more specificially:
%                  the e^-2 width of the intensity of the Gaussian pulse  =  4*std_i  =  FWHM_i*sqrt(2/log(2)),
%                  corresponding to the e^-1 width of the pulse field  =  sqrt(8)*std_f  =  FWHM_f/sqrt(log(2))
%     centerWavelength: the central wavelength of the wave in vaccuum, default: 1 micrometer
%     propagationSpeed: the propagation speed of the pulse in the medium
%
% spectralWidthMeasuredInWavelength: the e^-2 width of the intensity of the Gaussian spectrum in wavelength [m]
% widthOfIntensityInHz: the e^-2 width of the intensity of the Gaussian spectrum in Hz
%
% See also: calcGaussianSpectrumWidthInWavelengthFromPulseLength, and calcGaussianPulseLengthFromWidthInWavelength
function [spectralWidthMeasuredInWavelength, spectralWidthMeasuredInHz] = calcGaussianSpectrumWidthInWavelengthFromPulseDuration(pulseDurationMeasured, centerWavelength, propagationSpeed)
    if nargin<1 || isempty(pulseDurationMeasured)
        pulseDurationMeasured = 33*Const.u.femto; % 33 fs
    end
    if nargin<2 || isempty(centerWavelength)
        centerWavelength = Const.u.micro;
    end
    if nargin<3 || isempty(propagationSpeed)
        propagationSpeed = Const.c; % [a.u.]; Const.c, but cancels out anyway
    end
    
    centerFrequency = propagationSpeed./centerWavelength; % [Hz] ~ propagationSpeed
    
    % convert from the time to the frequency domain
    spectralWidthMeasuredInHz = calcGaussianSpectrumWidthInHzFromPulseDuration(pulseDurationMeasured);% [Hz] ~ propagationSpeed
    
    % Assume we only consider non-negative frequencies
    spectralIntervalMeasuredInHz = max(0,centerFrequency+spectralWidthMeasuredInHz*[-0.5 0.5]); % ~ propagationSpeed
    % convert frequencies to wavelengths
    spectralIntervalMeasuredInWavelength = propagationSpeed./spectralIntervalMeasuredInHz(end:-1:1); % ~ 1
    
    % finally, determine the width
    spectralWidthMeasuredInWavelength = diff(spectralIntervalMeasuredInWavelength); % ~ 1
    
%     % or simplified:
%     spectralWidthMeasuredInHzNormalized = (stdPerIntensityMeasurementInterval^2/(4*pi))*centerWavelength./pulseLengthMeasured;
%     spectralIntervalMeasuredInWavelengthSimplified = centerWavelength.*spectralWidthMeasuredInHzNormalized./max(0,1-(spectralWidthMeasuredInHzNormalized/2).^2);
%     spectralIntervalMeasuredInWavelengthApprox = (stdPerIntensityMeasurementInterval^2/(4*pi))*centerWavelength.^2./pulseLengthMeasured;

end
