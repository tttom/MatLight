% [pulseLengthMeasured pulseDurationMeasured] = calcGaussianPulseLengthFromSpectrumWidthInWavelength(spectralWidthMeasuredInWavelength, centerWavelength, stdPerIntensityMeasurementInterval)
%
% Determines the duration and length of Gaussian pulses with a given spectral width.
%
% Inputs:
%     spectralWidthMeasuredInHz: the width of the Gaussian spectrum in Hz, more specificially:
%                  the e^-2 width of the intensity of the Gaussian pulse  =  4*std_i  =  FWHM_i*sqrt(2/log(2)),
%                  corresponding to the e^-1 width of the pulse field  =  sqrt(8)*std_f  =  FWHM_f/sqrt(log(2))
%     centerWavelength: the central wavelength of the wave in vaccuum, default: 1 micrometer
%     stdPerIntensityMeasurementInterval: the number of standard deviation
%                   that fit within the two-sided measurement interval. 2 for e^-1 measurement, 4 for e^-2 (default)
%
% pulseLengthMeasured: the e^-2 length of the intensity of the Gaussian pulse in seconds
% pulseDurationMeasured: the e^-2 duration of the intensity of the Gaussian pulse in seconds
%
% See also: calcGaussianSpectrumWidthInWavelengthFromPulseLength, and calcGaussianPulseLengthFromWidthInWavelength
function [pulseLengthMeasured, pulseDurationMeasured] = calcGaussianPulseLengthFromSpectrumWidthInWavelength(spectralWidthMeasuredInWavelength, centerWavelength, stdPerIntensityMeasurementInterval)
    if nargin<1 || isempty(spectralWidthMeasuredInWavelength)
        spectralWidthMeasuredInWavelength = Const.u.nano; % [m]
    end
    if nargin<2 || isempty(centerWavelength)
        centerWavelength = Const.u.micro; % [m]
    end
    propagationSpeed = Const.c; % [m/s]
    if nargin<4 || isempty(stdPerIntensityMeasurementInterval)
        % How do we define the measured 'width' in terms of standard deviations?
        stdPerIntensityMeasurementInterval = 4; % e^-2 width of the intensity of the Gaussian pulse and Gaussian spectrum
    end
    
    spectralIntervalMeasuredInWavelength = centerWavelength + spectralWidthMeasuredInWavelength*[-0.5 0.5];
    spectralIntervalMeasuredInHz = propagationSpeed ./ spectralIntervalMeasuredInWavelength(end:-1:1);
    spectralWidthMeasuredInHz = diff(spectralIntervalMeasuredInHz);
    
    pulseDurationMeasured = calcGaussianPulseDurationFromSpectrumWidthInHz(spectralWidthMeasuredInHz, stdPerIntensityMeasurementInterval);
    
    pulseLengthMeasured = propagationSpeed * pulseDurationMeasured;
end
