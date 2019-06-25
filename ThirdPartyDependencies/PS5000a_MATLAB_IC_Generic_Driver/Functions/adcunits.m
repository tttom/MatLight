function adc_units = adcunits( time_units )
%ADCUNITS Returns a string representation of the time unit index.
%   adcunits(time_units) returns a string corresponding to the time unit index provided after
%   incrementing the index by 1.
%   
%   time_units can be one of the following:
%       
%       -1: ADC counts
%        0: Femtoseconds
%        1: Picoseconds
%        2: Nanoseconds
%        3: Microseconds
%        4: Milliseconds
%        5: Seconds
%
%   Applies only to the PicoScope 2000 Series devices using the ps2000.dll
%   driver.

    time_units = time_units + 1;
    %fprintf('time unit:  %d\n', time_units');
    
    switch(time_units)
        
        case 0
            
            adc_units = 'ADC';
            
        case 1
            
             adc_units = 'fs';
             
        case 2
            
             adc_units = 'ps';
            
        case 3
            
             adc_units = 'ns';
             
        case 4
            
             adc_units = 'us';
             
        case 5
            
             adc_units = 'ms';
             
        case 6
            
            adc_units = 's'
             
        otherwise
            
            adc_units = 'Not Known';

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    adcunits.m
%
% Copyright:   Pico Technology Limited 2013
%
% Author:      HSM
%
% Description:
%   This is a MATLAB script that returns a string representation of the 
%   time unit index.
%
% To call this function:
%
%   Type 'adc_units(time_units)' where 'time_units' is a value as per the 
%   case statement above.
%
%	Ensure that the location of this file is in your MATLAB Path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

