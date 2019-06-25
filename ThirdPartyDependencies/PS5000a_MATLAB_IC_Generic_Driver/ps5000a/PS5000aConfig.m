%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    PS5000aConfig.m
%
% Copyright:   Pico Technology Limited 2013
%
% Author:      HSM
%
% Description:
%   
%   Contains configuration data for setting parameters for a PicoScope 5000
%   Series Oscilloscope device.
%
%   Run this script in the MATLAB environment prior to connecting to the 
%   device.
%
%   This file can be edited to suit the application requirements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETUP PATH

%addpath('C:\Pico SDK\PS5000asdk_rx_x_x_x\') %Edit to specify location of dll files if required.
addpath('..');              % PicoStatus.m & PicoConstants.m
addpath('..\Functions');    % Common functions

%% LOAD ENUMS AND STRUCTURES

% Load in enumerations and structure information - DO NOT EDIT THIS SECTION
[ps5000aMethodinfo, ps5000aStructs, ps5000aEnuminfo, ps5000aThunkLibName] = ps5000aMFile; 

%% PICOSCOPE SETTINGS
% Define Settings for PicoScope 5000 series - if using tmtool, these can be
% accessed from the MATLAB environment by calling:
% evalin('base', 'variable_name');

