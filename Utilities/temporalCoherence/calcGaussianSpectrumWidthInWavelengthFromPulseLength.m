% [spectralWidthMeasuredInWavelength spectralWidthMeasuredInHz] = calcGaussianSpectrumWidthInWavelengthFromPulseLength(pulseLengthMeasured, centerWavelength)
%
% Determines the bandwidth for Gaussian pulses of a certain length.
%
% Inputs:
%     pulseLengthMeasured: the length of the pulse in the medium, more specificially:
%                  the e^-2 width of the intensity of the Gaussian pulse  =  4*std_i  =  FWHM_i*sqrt(2/log(2)),
%                  corresponding to the e^-1 width of the pulse field  =  sqrt(8)*std_f  =  FWHM_f/sqrt(log(2))
%     centerWavelength: the central wavelength of the wave in vaccuum, default: 1 micrometer
%
% spectralWidthMeasuredInWavelength: the e^-2 width of the intensity of the Gaussian spectrum in wavelength [m]
% widthOfIntensityInHz: the e^-2 width of the intensity of the Gaussian spectrum in Hz
%
% See also: calcGaussianSpectrumWidthInWavelengthFromPulseDuration, and calcGaussianPulseLengthFromWidthInWavelength
function [spectralWidthMeasuredInWavelength, spectralWidthMeasuredInHz] = calcGaussianSpectrumWidthInWavelengthFromPulseLength(pulseLengthMeasured, centerWavelength)
    if nargin<1 || isempty(pulseLengthMeasured)
        pulseLengthMeasured = 10e-6; % 1 meter
    end
    if nargin<2 || isempty(centerWavelength)
        centerWavelength = Const.u.micro;
    end
    propagationSpeed = Const.c; % cancels out
      
    pulseDurationMeasured = pulseLengthMeasured; % [s] % ~ 1./propagationSpeed, removed, to be added later again
    
    [spectralWidthMeasuredInWavelength, spectralWidthMeasuredInHz] = calcGaussianSpectrumWidthInWavelengthFromPulseDuration(pulseDurationMeasured, centerWavelength, 1.0);
    
    spectralWidthMeasuredInHz = propagationSpeed * spectralWidthMeasuredInHz; % [Hz], added propagationSpeed factor again
    
%     % or simplified:
%     stdPerIntensityMeasurementInterval = 4;
%     spectralWidthMeasuredInHzNormalized = (stdPerIntensityMeasurementInterval^2/(4*pi))*centerWavelength./pulseLengthMeasured;
%     spectralIntervalMeasuredInWavelengthSimplified = centerWavelength.*spectralWidthMeasuredInHzNormalized./max(0,1-(spectralWidthMeasuredInHzNormalized/2).^2);
%     spectralIntervalMeasuredInWavelengthApprox = (stdPerIntensityMeasurementInterval^2/(4*pi))*centerWavelength.^2./pulseLengthMeasured;

end
