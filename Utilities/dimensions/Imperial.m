% 
% A class representing the imperial units
%
classdef Imperial < English
    properties (Constant)
        mile=1609.344; % m
        gallon=4.54609*Unit.liter; % m^3
        fluidOunce=Imperial.gallon/160; % m^3
        pint=20*Imperial.fluidOunce; % m^3
    end
end
