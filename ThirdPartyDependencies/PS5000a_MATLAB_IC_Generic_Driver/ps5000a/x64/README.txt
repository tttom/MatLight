PicoScope 5000 Series MATLAB Generic Instrument Driver Release Notes
====================================================================

Driver Version: 1.1.9 Beta Release, 18th October 2013

License
-------

The driver and supporting files in this zip file are subject to the License terms in the PicoScope 5000 Series Programmer's Guide which can be found in the User Manuals section of our Support pages and in this zip file:

http://www.picotech.com/document/pdf/ps5000abpg.en-1.pdf

Installation
------------

For the underlying dlls required please contact support@picotech.com stating whether you have a 32 or 64-bit version of MATLAB.

Please add the location of the dlls to your MATLAB path (use the 'addpath' command if required).  

The mex -setup command may need to be run on your PC in order to load the dll files.

Ensure that the ps5000a.dll is named in lower case or rename it.

If PicoScope 6 has not been installed on your PC, please ensure it is for the USB driver.

The root folder of this zip file contains the following folders and files:

Functions 	- helper functions to convert data
PS5000a    	- Containing the instrument control driver (picotech_ps5000a_generic.mdd), information extracted from the header files, PicoScope 5000 specific constants, a configuration file & scripts.

PicoConstants.m - Constant values used for Pico devices - DO NOT EDIT
PicoStatus.m    - Status code definitions - DO NOT EDIT

The driver was created using MATLAB 2012b (32-bit for Windows) with Instrument Control Box v3.2.

The PS5000aConfig.m will setup relative paths to the Functions and root directory when run - ensure that this script is called before creating an icdevice object.


64-bit drivers
--------------

64-bit dlls are provided but have not been tested.

Copy the ps5000aMFile and ps5000aWrapMFile from the ps5000a/x64 directory to the ps5000a directory.


Supported Devices
-----------------

PicoScope 5242A/B & 5442A/B
PicoScope 5243A/B & 5443A/B
PicoScope 5244A/B & 5444A/B


Functions & Properties
----------------------

The Instrument Control driver contains a number of functions and properties which can be used to configure the device and collect data. 

For functions beginning with 'ps5000a' e.g. ps5000aSetChannel, please refer to the corresponding underlying function in the main Programmer's guide for further information on parameters and function descriptions. 

New functions that simplify calls to the underlying driver e.g. setSigGenBuiltInSimple have also been included in the driver.

The resetDevice() function can be used to set the driver back to it's original state which includes Channel A and B (as well as C and D for 4-channel devices) being set to a range of 5V.

Help text for some of the functions and properties can be found by loading the driver by using the 'instrhelp' command or in the Test and Measurement tool (tmtool) using the .


Examples Provided
-----------------

PS5000a_IC_Generic_Driver_Block.m 		- Captures a block of data using the driver's default settings and a simple trigger.

PS5000a_IC_Generic_Driver_Block_FFT.m 		- Captures a block of data using the driver's default settings and a simple trigger, then plots an FFT of the signal using an example from the MathWorks website.

PS5000a_IC_Generic_Driver_Rapid_Block.m  	- Captures a rapid block of data using a simple trigger.

PS5000a_IC_Generic_Driver_Rapid_Block_Plot3D.m  - Captures a rapid block of data using a simple trigger and plots the captures in 3D.

PS5000a_IC_Generic_Driver_Streaming.m  		- Captures a streaming data using a simple trigger and plots the captures for 2 channels.

PS5000a_Generic_Driver_Sig_Gen.m 		- Shows different examples of using the in-built function generator with specific waveforms.

These scripts can be run from the MATLAB command window or editor.


Test and Measurement Tool
-------------------------

The driver can be used with the Test and Measurement Tool by typing 'tmtool' at the MATLAB command prompt. 

Ensure that the current directory is set to PS5000a and that the path to the 'Functions' and root directory is set as described above. The PS5000aConfig.m file MUST be run PRIOR to loading the driver.

Create a new device object, and load in the driver file. Click 'Connect' to make the connection to the scope and it should then be ready for data collection.

Variables from the MATLAB workspace can be passed as arguments using the following command:

evalin('base', 'variable_name')

where the variable name must be in quotes.

Ensure that the device is disconnected when finished to avoid a driver lockup.


Unsupported Features
--------------------

The following features are not tested fully in this release of the driver:

- Advanced Triggering

Fixes
-----

* Minor fixes in functions
* Improved buffer setup for getRapidBlockData function

Known Issues
------------

* Depending on the PC, plotting data while collecting streaming data may slow down data collection.
* Setting up data buffers for rapid block data may take some time depending on the PC.
 

Contact
-------

Please send any questions, comments and feedback to support@picotech.com