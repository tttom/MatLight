% 
% A class representing the old english units
%
classdef English
    properties (Constant)
        % Imperial / American units
        inch=25.4e-3; % m
        feet=12*English.inch; % m
        yard=36*English.inch; % m
        nauticalMile=1852; %m
        pound=0.45359237; % kg
        ounce=31.1034768e-3; % kg
    end
end
