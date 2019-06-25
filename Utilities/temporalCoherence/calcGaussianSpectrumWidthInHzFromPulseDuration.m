% spectralWidthMeasuredInHz = calcGaussianSpectrumWidthInHzFromPulseDuration(pulseDurationMeasured,stdPerIntensityMeasurementInterval)
%
% Determines the bandwidth for Gaussian pulses of a certain duration.
%
% Inputs:
%     pulseDurationMeasured: the length of the pulse in time, more specificially:
%                  the e^-2 width of the intensity of the Gaussian pulse  =  4*std_i  =  FWHM_i*sqrt(2/log(2)),
%                  corresponding to the e^-1 width of the pulse field  =  sqrt(8)*std_f  =  FWHM_f/sqrt(log(2))
%     stdPerIntensityMeasurementInterval: the number of standard deviation
%     that fit within the two-sided measurement interval. 2 for e^-1 measurement, 4 for e^-2 (default)
%
% widthOfIntensityInHz: the e^-2 width of the intensity of the Gaussian spectrum in Hz
%
% See also: calcGaussianSpectrumWidthInWavelengthFromPulseLength, and calcGaussianPulseLengthFromWidthInWavelength
function spectralWidthMeasuredInHz = calcGaussianSpectrumWidthInHzFromPulseDuration(pulseDurationMeasured,stdPerIntensityMeasurementInterval)
    if nargin<1 || isempty(pulseDurationMeasured),
        pulseDurationMeasured = 33*Const.u.femto; % 33 fs
    end
    if nargin<2 || isempty(stdPerIntensityMeasurementInterval),
        % How do we define the measured 'width' in terms of standard deviations?
        stdPerIntensityMeasurementInterval = 4; % e^-2 width of the intensity of the Gaussian pulse and Gaussian spectrum
    end
    
    % internally we'll work with fields instead of intensities
    fieldStdPerMeasurementInterval = stdPerIntensityMeasurementInterval/sqrt(2);
    
    % convert from the time to the frequency domain
    pulseDurationFieldStd = pulseDurationMeasured./fieldStdPerMeasurementInterval; % [s]
    spectralWidthFieldStd = (1/(2*pi))./pulseDurationFieldStd; % [Hz] ~ propagationSpeed
    spectralWidthMeasuredInHz = fieldStdPerMeasurementInterval*spectralWidthFieldStd; % [Hz] ~ propagationSpeed
    
end
