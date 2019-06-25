% NKT Photonics SuperK class
%
%
%  Example usage
% 
%     % Initialize
%     source=SuperK();         % Create the object 'source' of the class SuperK
%     source.targetPower=0.50; % Set the power to 50%
% 
%     % Operate like this
%     source.setWavelengths(633e-9);  % He-Ne on full blast
%     source.setWavelengths(532e-9,0.5);  % freq-doubled Nd-YAG on half amplitude modulation
%     source.setWavelengths([633 532]*1e-9,[1 0.5]); % Turn on both at once
%     source.setWavelengths([633 532]*1e-9,0.5); % Set both to 50% amplitude
%     source.addWavelengths([450e-9 600e-9],1); % Add two extra wavelengths on full power.
% 
%     source % Get some readings
%     [wavelengths relativePowers gains]=source.getWavelengths()  % Get the settings back if you forgot
% 
%     % Shut down and free connection
%     source.delete();  % Shut everything down

classdef SuperK < LightSource
    % NKT Photonics SuperK class
    %
    properties(Access = private)
        serialPort;
        computerAddress;
        extremeSystemAddress;
        RFAddress;
        selectAddress;
        hasAOTF;
        intensityCorrectionFunctor=@(wavelengths) ones(size(wavelengths));
    end
    properties(Dependent = true)
        interlock;
        emission;
        constantPower;
        targetPower;
        RFPowerOn;
        intensityCorrection;
    end
    properties(SetAccess=private)
        maxNumberOfWavelengths=8;
        minWavelength; % [m]
        maxWavelength; % [m]
        temperature;   % [K]
        temperatureRF;   % [K]
        keySwitchOn; % 0: off, 1:on
        
        wavelengths; % [m]
        relativePowers; % fraction
        gains; % fraction
    end
    
    methods
        function source=SuperK()
            source.computerAddress=hex2dec('66');
            source.extremeSystemAddress=hex2dec('0F');
            source.RFAddress=hex2dec('06'); % hex2dec('10')+switchSetting;
            source.selectAddress=hex2dec('10');
            
            serialInfo = instrhwinfo('serial');
            potentialPorts=serialInfo.SerialPorts; %AvailableSerialPorts;
            moduleType=[];
            for portIdx=1:length(potentialPorts)
                port=potentialPorts{portIdx};
                openSerialPortObjects=instrfindall('Port',port,'Status','open');
                for (openSerialPortObject=openSerialPortObjects)
                    fclose(openSerialPortObject);
                    delete(openSerialPortObject);
                end
                source.serialPort=serial(port,'BaudRate',115200,'DataBits',8,'StopBits',1,'Parity','none');
                fopen(source.serialPort);
                try
                    moduleType=source.readRegister(source.extremeSystemAddress,hex2dec('61'),.2);
                    if (moduleType==hex2dec('60'))
                        logMessage('Talking to SuperK over %s.',port);
                        break;
                    end
                catch Exc
                    logMessage(['DEBUG: ''' Exc.message ''' on line ' num2str(Exc.stack(1).line) ' in file ' Exc.stack(1).file]);
                end
                fclose(source.serialPort);
                source.serialPort=[];
            end
            if (moduleType==hex2dec('60'))
                source.initialize();
            else
                logMessage('Error: cannot connect to SuperK over serial %d ports.',length(potentialPorts));
            end
        end
        function set.interlock(source,value)
            source.writeRegister(source.extremeSystemAddress,hex2dec('32'),~value); % Interlock reset / set
        end
        function value=get.interlock(source)
            [ign value]=source.getStatus();
        end
        function value=get.keySwitchOn(source)
            [ign ign ign interlockLoopOff]=source.getStatus();
            value=~interlockLoopOff;
        end
        function set.emission(source,value)
            source.interlock=false; % Interlock reset
            source.writeRegister(source.extremeSystemAddress,hex2dec('30'),value*3); % Emission off / on
            source.constantPower=true; % Constant Power
        end
        function value=get.emission(source)
            response=source.readRegister(source.extremeSystemAddress,hex2dec('30'));
            value=response>0;
        end
        function set.RFPowerOn(source,value)
            if source.hasAOTF,
                source.writeRegister(source.RFAddress,hex2dec('30'),value); % RF Power off / on
                source.writeRegister(source.RFAddress,hex2dec('31'),2^1+2^0); % Setup bits: Optimal power table, and filter temperature compensation
                source.minWavelength=double(source.readRegister(source.RFAddress,hex2dec('34')))*1e-12;
                source.maxWavelength=double(source.readRegister(source.RFAddress,hex2dec('35')))*1e-12;
            else
                logMessage('No AOTF found at initialization.');
            end
        end
        function value=get.RFPowerOn(source)
            if source.hasAOTF,
                response=source.readRegister(source.RFAddress,hex2dec('30'));
                value=response>0;
            else
                logMessage('No AOTF found at initialization.');
            end
        end
        function set.constantPower(source,value)
            source.writeRegister(source.extremeSystemAddress,hex2dec('31'),value); % Constant Power off / on
        end
        function value=get.constantPower(source)
            response=source.readRegister(source.extremeSystemAddress,hex2dec('31'));
            value=response>0;
        end
        function value=get.temperature(source)
            value=double(source.readRegister(source.extremeSystemAddress,hex2dec('11')))/10+273.15;
        end
        function value=get.temperatureRF(source)
            if source.hasAOTF,
                value=double(source.readRegister(source.RFAddress,hex2dec('38')))/10+273.15;
            else
                logMessage('No AOTF found at initialization.');
            end
        end
        function set.targetPower(source,newPowerFraction)
            source.interlock = false; % Interlock reset
            source.emission = (newPowerFraction>0); % Emission ON when power > 0
            source.constantPower = true; % Constant Power
            source.writeRegister(source.extremeSystemAddress,hex2dec('37'),uint32(newPowerFraction*1000)); % Set power in permille
        end
        function targetPower=get.targetPower(source)
            response=source.readRegister(source.extremeSystemAddress,hex2dec('37'));
            targetPower=double(response)/1000.0;
        end
        function setWavelengths(source,wavelengths,relativePowers,gains)
            if (nargin<3)
                relativePowers=[];
            end
            if (length(relativePowers)<length(wavelengths))
                if (isempty(relativePowers))
                    relativePowers=1; % Maximum power
                end
                relativePowers(end+1:length(wavelengths))=relativePowers(end); % Set all undefined powers to the last one defined
            end
            if (nargin<4)
                gains=[];
            end
            if (length(gains)<length(wavelengths))
                if (isempty(gains))
                    gains=1; % Maximum gain
                end
                gains(end+1:length(wavelengths))=gains(end); % Set all undefined gains to the last one defined
            end
            % Adjust the intensity so the measured intensity is the requested
            relativePowers=relativePowers.*source.intensityCorrectionFunctor(wavelengths);
            % Make sure no negative powers or gains are presented
            relativePowers=max(0,relativePowers);
            gains=max(0,gains);
            % Zero pad
            relativePowers(source.maxNumberOfWavelengths+1)=0;
            % One pad
            gains(source.maxNumberOfWavelengths+1)=1;
               
            if source.hasAOTF,
                source.RFPowerOn=any(relativePowers>0);
                for (wavelengthIdx=1:source.maxNumberOfWavelengths)
                    if (wavelengthIdx<=length(wavelengths))
                        if (wavelengths(wavelengthIdx)~=source.wavelengths(wavelengthIdx))
                            source.writeRegister(source.RFAddress,hex2dec('90')+wavelengthIdx-1,wavelengths(wavelengthIdx)*1e12); % Set wavelength
                            source.wavelengths(wavelengthIdx)=wavelengths(wavelengthIdx);
                        end
                    end
                    if (relativePowers(wavelengthIdx)~=source.relativePowers(wavelengthIdx))
                        source.writeRegister(source.RFAddress,hex2dec('B0')+wavelengthIdx-1,relativePowers(wavelengthIdx)*1000); % Set amplitude
                        source.relativePowers(wavelengthIdx)=relativePowers(wavelengthIdx);
                    end
                    if (gains(wavelengthIdx)~=source.gains(wavelengthIdx))
                        source.writeRegister(source.RFAddress,hex2dec('C0')+wavelengthIdx-1,gains(wavelengthIdx)*1000); % Set gain
                        source.gains(wavelengthIdx)=gains(wavelengthIdx);
                    end
                end
            else
                logMessage('No AOTF found at initialization.');
            end
        end
        function [wavelengths, relativePowers, gains]=getWavelengths(source)
            wavelengths=source.wavelengths;
            relativePowers=source.relativePowers;
            gains=source.gains;
            
            % Return only the active ones
            idx=find(relativePowers>0 & gains>0,true,'last');
            wavelengths=wavelengths(1:idx);
            relativePowers=relativePowers(1:idx);
            % Return the set powers according to the new calibration.
            relativePowers=relativePowers./source.intensityCorrectionFunctor(wavelengths);
            gains=gains(1:idx);
        end
        function addWavelengths(source,extraWavelengths,extraPowers,extraGains)
            if (nargin<3 || isempty(extraPowers))
                extraPowers=1;
            end
            if (length(extraPowers)<length(extraWavelengths))
                extraPowers(end+1:length(extraWavelengths))=extraPowers(end); % Replicate the last value for all new wavelengths
            end
            if (nargin<4 || isempty(extraGains))
                extraGains=1;
            end
            if (length(extraGains)<length(extraWavelengths))
                extraGains(end+1:length(extraWavelengths))=extraGains(end); % Replicate the last value for all new wavelengths
            end
            % Append the new wavelengths
            [currentWavelengths currentRelativePowers currentGains]=source.getWavelengths();
            newWavelengths=[currentWavelengths extraWavelengths];
            newRelativePowers=[currentRelativePowers extraPowers];
            newGains=[currentGains extraGains];
            % Check the length
            if (length(newWavelengths)>source.maxNumberOfWavelengths)
                logMessage('Attempting to set %d wavelengths, restricting it to %d.',[length(newWavelengths) source.maxNumberOfWavelengths]);
                newWavelengths=newWavelengths(1:source.maxNumberOfWavelengths);
                newRelativePowers=newRelativePowers(1:source.maxNumberOfWavelengths);
                newGains=newGains(1:source.maxNumberOfWavelengths);
            end
            %Update the wavelengths
            source.setWavelengths(newWavelengths,newRelativePowers,newGains);
        end
        function setIntensityCorrection(source,wavelengths,intensityCorrectionValues)
            if (nargin<3)
                if (isa(wavelengths,'function_handle'))
                    source.intensityCorrection=wavelengths;
                else
                    logMessage('Invalid argument specified to setIntensityCorrection(source,wavelengths,measuredValues). Either specify a set of wavelengths followed by the measured values, or specify a function that returns the correction for each wavelength.');
                end
            else
                source.intensityCorrection=@(w) interp1(wavelengths,intensityCorrectionValues,w,'cubic',1);
%                 intensityCorrectionSmoothingSpline=spaps(x.',y.',100*1e-9);
%                 source.intensityCorrection=@(w) fnval(intensityCorrectionSmoothingSpline,w);
            end
        end
        function set.intensityCorrection(source,intensityCorrectionFunctor)
            source.intensityCorrectionFunctor=intensityCorrectionFunctor;
        end
        function intensityCorrectionFunctor=get.intensityCorrection(source)
            intensityCorrectionFunctor=source.intensityCorrectionFunctor;
        end
        function delete(source)
            source.deinitialize();
            
            fclose(source.serialPort);
            clear source.serialPort;
        end
        
        function hexString=convertToHex(source,message)
            hexString=dec2hex(message(:).').';
            hexString=hexString(:).';
        end
    end
    methods(Access=private)
        function initialize(source)
            % Check if AOTF present
            try
                moduleType=source.readRegister(source.RFAddress,hex2dec('61'),.2);
                if (moduleType==hex2dec('60'))
                    logMessage('RF generator found, AOTF installed.');
                    source.hasAOTF = true;
                else
                    logMessage('No RF generator found, AOTF not installed.');
                    source.hasAOTF = false;
                end
            catch Exc
                logMessage(['DEBUG: ''' Exc.message ''' on line ' num2str(Exc.stack(1).line) ' in file ' Exc.stack(1).file]);
            end
                
            % Initial settings
            source.interlock=false;
            if (source.hasAOTF),
                loadWavelengths(source);
            end
            if (~source.keySwitchOn)
                logMessage('Note: The keyswitch on the front panel is not turned to ON yet!');
            end
        end
        function deinitialize(source)
            source.targetPower=0.0;
            source.emission=false; % Emission OFF
            if source.hasAOTF,
                source.setWavelengths([]);
                source.RFPowerOn=false;
            end
        end
        function [wavelengths, relativePowers, gains]=loadWavelengths(source)
            wavelengths=[];
            relativePowers=[];
            gains=[];
            if source.hasAOTF,
                for (wavelengthIdx=1:source.maxNumberOfWavelengths)
                    nextWavelength=source.readRegister(source.RFAddress,hex2dec('90')+wavelengthIdx-1);
                    wavelengths(end+1)=nextWavelength(1);
                    nextPower=source.readRegister(source.RFAddress,hex2dec('B0')+wavelengthIdx-1);
                    relativePowers(end+1)=nextPower(1);
                    nextGain=source.readRegister(source.RFAddress,hex2dec('B0')+wavelengthIdx-1);
                    gains(end+1)=nextGain(1);
                end
            end
            wavelengths=double(wavelengths)*1e-12;
            relativePowers=double(relativePowers)/1000;
            
            source.wavelengths=wavelengths;
            source.relativePowers=relativePowers;
            source.gains=gains;
        end
        function [emission, interlock, interlockPowerFailure, interlockLoopOff, externalDisable, supplyVoltageLow, moduleTempRange, USBLogErrorPresent, errorCodePresent]=getStatus(source)
            response=source.readRegister(source.extremeSystemAddress,hex2dec('66'));
            emission=bitand(response,2^0)>0;
            interlock=bitand(response,2^1)==0;
            interlockPowerFailure=bitand(response,2^2)>0;
            interlockLoopOff=bitand(response,2^3)>0;
            externalDisable=bitand(response,2^4)>0;
            supplyVoltageLow=bitand(response,2^5)>0;
            moduleTempRange=bitand(response,2^6)>0;
            USBLogErrorPresent=bitand(response,2^14)>0;
            errorCodePresent=bitand(response,2^15)>0;
        end
        %% Register access
        function writeRegister(source,destinationModule,register,value,timeOut)
            if (nargin<5)
                timeOut=10;
            end
            message=zeros(1,8,'uint8');
            message(1)=destinationModule;
            message(2)=source.computerAddress; % Reply-To port
            message(3)=hex2dec('05'); % Write
            message(4)=register;
            message(5:8)=typecast(uint32(value),'uint8');
            messageWithCRC=source.addCRCToMessage(message);
            escapedAndFramedMessage=source.escapeAndFrameMessage(messageWithCRC);
%             source.convertToHex(escapedAndFramedMessage)
            fwrite(source.serialPort,escapedAndFramedMessage,'uint8');
                        
            % Read ACK
            escapedAndFramedMessage=source.readTelegram(timeOut);
            if (~isempty(escapedAndFramedMessage))
%                 source.convertToHex(escapedAndFramedMessage)
                messageWithCRC=source.unFrameAndUnEscapeMessage(escapedAndFramedMessage);
                [message crcError]=source.checkAndRemoveCRCFromMessage(messageWithCRC);
                if (crcError)
                    logMessage('Warning: CRC error after receiving message.');
                end
                if (message(1)~=source.computerAddress)
                    logMessage('Warning: received message for %d instead of %d.',[message(1) source.computerAddress]);
                end
                if (message(2)~=destinationModule)
                    logMessage('Warning: received message from %d instead of %d.',[message(1) destinationModule]);
                end
                if (message(3)~=3)
                    logMessage('Warning: received message of type %d instead of 3.',message(3));
                end
            else
                logMessage('Did not receive acknowlegement response.');
            end
        end
        function value=readRegister(source,destinationModule,register,timeOut)
            if (nargin<4)
                timeOut=10;
            end
            % Send READ request
            message=zeros(1,5,'uint8');
            message(1)=destinationModule;
            message(2)=source.computerAddress; % Reply-To port
            message(3)=hex2dec('04'); % Read
            message(4)=register;
            messageWithCRC=source.addCRCToMessage(message);
            escapedAndFramedMessage=source.escapeAndFrameMessage(messageWithCRC);
%             logMessage('Request');
%             source.convertToHex(escapedAndFramedMessage)
            fwrite(source.serialPort,escapedAndFramedMessage,'uint8');
            
            % Read response
            escapedAndFramedMessage=source.readTelegram(timeOut);
            if (~isempty(escapedAndFramedMessage))
%                 logMessage('Response');
%                 source.convertToHex(escapedAndFramedMessage)
                messageWithCRC=source.unFrameAndUnEscapeMessage(escapedAndFramedMessage);
                [message crcError]=source.checkAndRemoveCRCFromMessage(messageWithCRC);
                if (crcError)
                    logMessage('Warning: CRC error after receiving message.');
                end
                if (message(1)~=source.computerAddress)
                    logMessage('Warning: received message for %d instead of %d.',[message(1) source.computerAddress]);
                end
                if (message(2)~=destinationModule)
                    logMessage('Warning: received message from %d instead of %d.',[message(1) destinationModule]);
                end
                if (message(3)~=8)
                    logMessage('Warning: received message of type %d instead of 8.',message(3));
                end
                payload=message(5:end);

                % Read
                switch length(payload)
                    case 0
                        value=[];
                    case 1
                        value=payload;
                    case 2
                        value=typecast(payload,'uint16');
                    case 3
                        value=typecast([0 payload],'uint32');
                        logMessage('Warning: byte missing.');
                    case 4
                        value=typecast(payload,'uint32');
                    otherwise
                        value=typecast(payload(1:(4*floor(end/4))),'uint32');
                        if (mod(length(payload),4)>0)
                            logMessage('Warning: dropping %d bytes.',mod(length(payload),4));
                        end
                end
            else
                value=[];
            end
        end
        %% Communication
        function telegram=readTelegram(source,timeOut)
            if (nargin<2)
                timeOut=10;
            end
            telegram=uint8([]);
            startTime=clock();
            while (isempty(timeOut) || etime(clock(),startTime)<timeOut) && (length(telegram)<1 || telegram(end)~=hex2dec('0A'))
                if (source.serialPort.BytesAvailable>0)
                    newByte=fread(source.serialPort,1,'uint8');
                    if (length(telegram)>0 || newByte==hex2dec('0D'))
                        telegram(end+1)=newByte;
                    else
                        logMessage('Warning: unknown bytes received before telegram!');
                    end
                else
                    pause(0.010); % Take a little break before trying again
                end
            end
        end
        
        function escapedAndFramedMessage=escapeAndFrameMessage(source,message)
            escapedAndFramedMessage=zeros(1,1,'uint8');
            escapedAndFramedMessage=hex2dec('0D');
            for (idx=1:length(message))
                nextByte=message(idx);
                switch(nextByte)
                    case hex2dec('0A')
                        escapedAndFramedMessage(end+1)=hex2dec('5E');
                        escapedAndFramedMessage(end+1)=hex2dec('4A');
                    case hex2dec('0D')
                        escapedAndFramedMessage(end+1)=hex2dec('5E');
                        escapedAndFramedMessage(end+1)=hex2dec('4D');
                    case hex2dec('5E')
                        escapedAndFramedMessage(end+1)=hex2dec('5E');
                        escapedAndFramedMessage(end+1)=hex2dec('9E');
                    otherwise
                        escapedAndFramedMessage(end+1)=nextByte;
                end
            end
            escapedAndFramedMessage(end+1)=hex2dec('0A');
        end
        function message=unFrameAndUnEscapeMessage(source,escapedAndFramedMessage)
            if (escapedAndFramedMessage(1)==hex2dec('0D'))
                escapedMessage=escapedAndFramedMessage(2:end);
            else
                escapedMessage=escapedAndFramedMessage;
                logMessage('Error: response did not not start with 0xOD.');
            end
            if (escapedMessage(end)==hex2dec('0A'))
                escapedMessage=escapedMessage(1:end-1);
            else
                logMessage('Error: missing data, response did not not end with 0xOA.');
            end

            message=zeros(1,0,'uint8');
            bytesToSkip=0;
            for (idx=1:length(escapedMessage))
                nextByte=escapedMessage(idx);
                if (bytesToSkip==0)
                    switch(nextByte)
                        case hex2dec('5E')
                            if (idx<length(escapedMessage))
                                switch (escapedMessage(idx+1)),
                                    case hex2dec('4A')
                                        message(end+1)=hex2dec('0A');
                                    case hex2dec('4D')
                                        message(end+1)=hex2dec('0D');
                                    case hex2dec('9E')
                                        message(end+1)=hex2dec('5E');
                                    otherwise
                                        logMessage('Error: unknown escape character %d',escapedMessage(idx+1));
                                end
                                bytesToSkip=1;
                            else
                                logMessage('Error: escape character found at end of message.');
                            end
                        otherwise
                            message(end+1)=nextByte;
                    end
                else
                    bytesToSkip=bytesToSkip-1;
                end
            end
        end

        function messageWithCRC=addCRCToMessage(source,message)
            messageWithCRC=[message source.crcCCITTxModem(message)];
        end

        function crc=crcCCITTxModem(source,messageBytes)
            % Update the CRC for transmitted and received data using
            % the CCITT 16bit algorithm (X^16 + X^12 + X^5 + 1).
            crc=uint16(0);
            for (byteIdx=1:length(messageBytes))
                crc = swapbytes(crc);
                crc = bitxor(crc,uint16(messageBytes(byteIdx)));
                crc = bitxor(crc,bitshift(bitand(crc,255),-4));
                crc = bitxor(crc,bitshift(bitand(crc,15),12));
                crc = bitxor(crc,bitshift(bitand(crc,255),5));
            end
            crc=typecast(swapbytes(crc),'uint8');
        end

        function [message, crcError]=checkAndRemoveCRCFromMessage(source,messageWithCRC)
            message=messageWithCRC(1:end-2);
            if (all(source.crcCCITTxModem(message)==typecast(messageWithCRC(end-1:end),'uint8')))
                crcError=false;
            else
                message=[]; % An error occured
                crcError=true;
            end
        end
    end
end
