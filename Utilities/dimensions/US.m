% 
% A class representing the American units
%
classdef US < English
    properties (Constant)
        mile=1609.347219; % m
        gallon=231*US.inch^3; % m^3
        fluidOunce=US.gallon/128; % m^3
        pint=16*US.fluidOunce; % m^3
    end
end
