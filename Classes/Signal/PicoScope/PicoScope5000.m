classdef PicoScope5000 < handle
    % PicoScope 5000 class
    %
    % constructor PicoScope5000()
    %
    %
    properties(Access=public,Dependent=true)
        channels;
        timeInterval;
    end
    properties(SetAccess=private,GetAccess=public,Dependent=true)
        maxChannels;
        maxSamples;
    end
    properties(Access=private)
        instrument;
        privateChannels;
        mFile;
    end
    methods(Access=private,Static=true)
        function prepareBase()
            % Litter base to please bad mdd file
            evalin('base','[ps5000aMethodinfo, ps5000aStructs, ps5000aEnuminfo, ps5000aThunkLibName] = ps5000aMFile();');
        end
        function cleanBase()
            % Clean up base
            evalin('base','clear ps5000aMethodinfo ps5000aStructs ps5000aEnuminfo ps5000aThunkLibName;');
        end
        function varargout=invokeWithBase(varargin)
            varargout=cell(nargout,1);
            PicoScope5000.prepareBase();
            [varargout{:}]=invoke(varargin{:});
            PicoScope5000.cleanBase();
        end
    end
    methods
        function ps=PicoScope5000()
            ps.mFile=struct();
            [ps.mFile.ps5000aMethodinfo, ps.mFile.ps5000aStructs, ps.mFile.ps5000aEnuminfo, ps.mFile.ps5000aThunkLibName] = ps5000aMFile();
            
            ps.prepareBase();
            ps.instrument=instrfind();
            if ~isempty(ps.instrument)
                ps.instrument=ps.instrument(strcmpi(ps.instrument.Type,'scope')&strcmpi(ps.instrument.Status,'open'));
            end
            %ps.instrument=instrfind('DriverName','picotech_ps5000a_generic.mdd','Status','open');
            if isempty(ps.instrument)
                ps.instrument=icdevice('picotech_ps5000a_generic.mdd');
                connect(ps.instrument);
            else
                ps.instrument=ps.instrument(1);
                if strcmpi(ps.instrument.Status,'open'),
                    %disconnect(ps.instrument);
                    invoke(ps.instrument, 'resetDevice');
                else
                    
                end
            end
            ps.cleanBase();
            
            ps.channels=[5 5];
            ps.timeInterval=1e-6;
            
            % Channel     : 0 (PS5000A_CHANNEL_A)
            % Threshold   : 1000 (mV)
            % Direction   : 2 (Rising)
            % Delay       : 0
            % Auto trigger: 0 (wait indefinitely)
            %status=ps.invokeWithBase(ps.instrument, 'setSimpleTrigger',0,1*1e3,2,0,0);
            status=ps.invokeWithBase(ps.instrument, 'setTriggerOff');
        end
        function nbChannels=get.maxChannels(ps)
            nbChannels=ps.instrument.channelCount;
        end
        function set.channels(ps,chs)
            % chs is a vector of four scalars indicating the DC when
            % positive, and the AC range when negative. When omitted or
            % zero, the channel is not enabled.
            allowedVoltageRanges=[ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_10MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_20MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_50MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_100MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_200MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_500MV,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_1V,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_2V,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_5V,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_10V,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_20V,...
                                  ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_50V];
            allowedVoltages=[0.010,0.020,0.050,0.100,0.200,0.500,1,2,5,10,20,50];
            for idx=1:ps.maxChannels,
                if idx<=numel(chs),
                    ch=chs(idx);
                else
                    ch=0;
                end
                coupling=ps.mFile.ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC;
                if ch~=0 && ~isnan(ch),
                    if ch<0,
                        ch=-ch;
                        coupling=ps.mFile.ps5000aEnuminfo.enPS5000ACoupling.PS5000A_AC;
                    end
                    [~, chI]=find(allowedVoltages>=ch,1,'first'); % Round up
                else
                    chI=[];
                end
                if ~isempty(chI)
                    status=ps.invokeWithBase(ps.instrument, 'ps5000aSetChannel', idx-1, PicoConstants.TRUE, coupling,...
                        allowedVoltageRanges(chI), 0.0);
                    ps.privateChannels(idx)=allowedVoltages(chI);
                    if coupling==ps.mFile.ps5000aEnuminfo.enPS5000ACoupling.PS5000A_AC,
                        ps.privateChannels(idx)=-ps.privateChannels(idx); % Indicate AC
                    end
                else
                    status=ps.invokeWithBase(ps.instrument, 'ps5000aSetChannel', idx-1, PicoConstants.FALSE, ps.mFile.ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC,...
                        ps.mFile.ps5000aEnuminfo.enPS5000ARange.PS5000A_5V, 0.0);
                    ps.privateChannels(idx)=0;
                end
            end
            
            % Update the time interval and sample resolution
            ps.timeInterval=ps.timeInterval;
        end
        function chs=get.channels(ps)
            chs=ps.privateChannels;
        end
        function maxSmpls=get.maxSamples(ps)
            [status,timeIntervalNanoSeconds,maxSmpls]=ps.invokeWithBase(ps.instrument,'ps5000aGetTimebase',ps.instrument.timebase,0);
        end
        function set.timeInterval(ps,dt)
            resolutions=[8 12 14 15 16]; %all potential bit rates
            nbChannels=sum(ps.channels~=0);
            switch (nbChannels)
                case 1
                    minInterval=[1e-9 2e-9 8e-9 8e-9 16e-9];
                case 2
                    minInterval=[2e-9 4e-9 8e-9 8e-9 Inf];
                otherwise
                    minInterval=[4e-9 8e-9 8e-9 Inf Inf];
            end
            validResolutions=(dt>=minInterval);
            maxBits=max(resolutions(validResolutions));
            if isempty(maxBits)
                maxBits=8;
            end
            
            % Update the device precision
            status=ps.invokeWithBase(ps.instrument, 'ps5000aSetDeviceResolution', maxBits);
            
            % Store the time base for the following operations
            if maxBits==12 || maxBits==16,
                dtThreshold=16e-9;
            else
                dtThreshold=8e-9;
            end
            if dt<dtThreshold,
                timeBase=max(0,round(log2(dt/1e-9)));
            else
                timeBase=round(dt/dtThreshold)+(log2(dtThreshold/1e-9)-1);
            end
            ps.instrument.timebase=timeBase; % Change the instrument settings too
                        
        end
        function timeInt=get.timeInterval(ps)
            [status,timeIntervalNanoSeconds]=invoke(ps.instrument,'ps5000aGetTimebase',ps.instrument.timebase,0);
            if status~=PicoStatus.PICO_OK,
                ps.instrument.timebase=65;
                [status,timeIntervalNanoSeconds]=invoke(ps.instrument,'ps5000aGetTimebase',ps.instrument.timebase,0);
            end

            timeInt=double(timeIntervalNanoSeconds)*1e-9;
        end
        function [data,times]=acquire(ps,nbSamples)
            set(ps.instrument, 'numPreTriggerSamples', 0);
            set(ps.instrument, 'numPostTriggerSamples', nbSamples);
            
            ps.prepareBase();
            status=invoke(ps.instrument,'runBlock',0);
            % Retrieve data values:
            %
            % start index       : 0
            % segment index     : 0
            % downsampling ratio: 1
            % downsampling mode : 0 (PS5000A_RATIO_MODE_NONE)
            [chA,chB,chC,chD,nbSamples,overflow]=invoke(ps.instrument,'getBlockData',0,0,1,0);
            chABCD={chA,chB,chC,chD}; clear chA chB chC chD;
            status=invoke(ps.instrument,'ps5000aStop');
            ps.cleanBase();
            % /maxADCValue
            
            % Prepare output
            activeChannels=find(ps.channels~=0.0);
            data=zeros(numel(activeChannels),nbSamples);
            for activeChannel=activeChannels,
                tmp=chABCD{activeChannel};
                data(activeChannel,:)=tmp(:).';
                clear tmp;
            end
            data=data*1e-3; % Convert to V
            if nargout>1,
                times=([1:double(nbSamples)]-1)*ps.timeInterval;
            end
        end
        function close(ps)
            ps.prepareBase();
            disconnect(ps.instrument);
            ps.cleanBase();
        end
        function delete(ps)
            ps.close();
        end
    end    
end
