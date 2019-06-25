function adc_count = mv2adc( mv, voltage_range, max_value )
% MV2ADC Converts milliVolt value to ADC count
%   mv2adc( mv, range, max_value ) returns an ADC count value corresponding
%   to the voltage range selected.
%
%   mv - the value in millivolts
%   voltage_range - the voltage range in millivolts
%   max_value - the ADC max value for the device

    adc_count = int32((int32(mv) * double(max_value)) / voltage_range);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    mv2adc.m
%
% Copyright:   Pico Technology Limited 2012
%
% Author:      HSM
%
% Description:
%   
%		This is a MATLAB script that converts millivolt values to ADC counts.
%
%       Ensure that the location of this file is in your MATLAB Path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

