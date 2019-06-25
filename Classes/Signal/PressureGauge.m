classdef PressureGauge < handle
    % Pressure Gauge class
    %
    % constructor PhaseSLM(serialPort)
    %
    %
    properties(SetAccess=private,GetAccess=public)
        serialPortNumber;
    end
    properties(Access=private)
        serialPort;
    end
    methods
        function pg=PressureGauge(serialPortNumber)
            if (nargin<1)
                serialPortNumber=[];
            end
            if isempty(serialPortNumber)
                availablePorts=getAvailableComPort();
                availablePorts=cellfun(@(e)str2double(e(4:end)),availablePorts).';
            end
            availableInstruments=instrfind('Type','serial');
            availableInstrumentsPort=availableInstruments.Port;
            if ~iscell(availableInstrumentsPort)
                availableInstrumentsPort={availableInstrumentsPort};
            end
            availableInstrumentsNumbers=cellfun(@(e) str2double(e(4:end)),availableInstrumentsPort);
            availableInstrumentsNumbers=availableInstrumentsNumbers(~isnan(availableInstrumentsNumbers)).';
            
            % Take the last one
            serialPortNumber=max([availablePorts availableInstrumentsNumbers]);
            
            portName=sprintf('COM%0.0f',serialPortNumber);
            % First close all
            availableInstrumentsWithThisPortNumber=instrfind('Type','serial','Port',portName);
            for instr=availableInstrumentsWithThisPortNumber,
                delete(instr);
            end
            % Open a new one
            pg.serialPort=serial(portName);
            
            set(pg.serialPort, 'BaudRate', 19200);
            set(pg.serialPort, 'FlowControl', 'hardware');
            set(pg.serialPort, 'Timeout', 2.0);
            fopen(pg.serialPort);
        end
        function spn=get.serialPortNumber(pg)
           spn=pg.serialPort.number;
        end
        function p=getPressure(pg)
            pressureString = query(pg.serialPort, 'IN_PV_1', '%s\r' ,'%s\n');
            p=str2double(pressureString(1:end-4));
            p=p*101325/760;
        end
        function close(pg)
            fclose(pg.serialPort);
        end
        function delete(pg)
            pg.close();
        end
    end    
end
