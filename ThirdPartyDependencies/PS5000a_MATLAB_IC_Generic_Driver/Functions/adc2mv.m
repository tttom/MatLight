function mv_value = adc2mv(adc_count, voltage_range, max_value)
%ADC2MV Converts raw ADC count value to millivolts
%   adc2mv(raw, voltageRange, maxValue) returns a millivolt
%   value corresponding to the ADC Count and the Voltage range set:
%
%       adc_count - the raw ADC value
%       voltage_range - the voltage range used for the channel (in millivolts)
%       max_value - the maximum ADC value for the device

    mv_value = (double(adc_count) * voltage_range) / max_value ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    adc2mv.m
%
% Copyright:   Pico Technology Limited 2012
%
% Author:      HSM
%
% Description:
%   This is a MATLAB script that converts ADC counts to millivolt values.
%
%	Ensure that the location of this file is in your MATLAB Path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%