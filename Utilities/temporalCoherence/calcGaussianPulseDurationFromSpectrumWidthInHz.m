% pulseDurationMeasured = calcGaussianPulseDurationFromSpectrumWidthInHz(spectralWidthMeasuredInHz,stdPerIntensityMeasurementInterval)
%
% Determines the pulse duration of Gaussian pulses for a given spectral width.
%
% Inputs:
%     spectralWidthMeasuredInHz: the width of the Gaussian spectrum in Hz, more specificially:
%                  the e^-2 width of the intensity of the Gaussian pulse  =  4*std_i  =  FWHM_i*sqrt(2/log(2)),
%                  corresponding to the e^-1 width of the pulse field  =  sqrt(8)*std_f  =  FWHM_f/sqrt(log(2))
%     stdPerIntensityMeasurementInterval: the number of standard deviation
%                   that fit within the two-sided measurement interval. 2 for e^-1 measurement, 4 for e^-2 (default)
%
% pulseDurationMeasured: the e^-2 width of the intensity of the Gaussian pulse in seconds
%
% See also: calcGaussianSpectrumWidthInWavelengthFromPulseLength, and calcGaussianPulseLengthFromWidthInWavelength
function pulseDurationMeasured = calcGaussianPulseDurationFromSpectrumWidthInHz(spectralWidthMeasuredInHz,stdPerIntensityMeasurementInterval)
    if nargin<1 || isempty(spectralWidthMeasuredInHz),
        spectralWidthMeasuredInHz = 1.2923e-007; % [Hz]
    end
    if nargin<2 || isempty(stdPerIntensityMeasurementInterval),
        % How do we define the measured 'width' in terms of standard deviations?
        stdPerIntensityMeasurementInterval = 4; % e^-2 width of the intensity of the Gaussian pulse and Gaussian spectrum
    end
    
    % internally we'll work with fields instead of intensities
    fieldStdPerMeasurementInterval = stdPerIntensityMeasurementInterval/sqrt(2);
    
    % convert from the frequency to the time domain
    spectralWidthFieldStd = spectralWidthMeasuredInHz./fieldStdPerMeasurementInterval; % [Hz] ~ propagationSpeed
    pulseDurationFieldStd = (1/(2*pi))./spectralWidthFieldStd; % [Hz] ~ propagationSpeed
    pulseDurationMeasured = fieldStdPerMeasurementInterval*pulseDurationFieldStd; % [s]
end
