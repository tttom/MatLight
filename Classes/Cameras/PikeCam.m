classdef PikeCam < Cam
    % Pike camera class
    %  
    properties (SetAccess = private)
        maxSize;
        minRegionOfInterestSize=[1 1]*2;
        bitsPerPixel;
    end
    properties (Access = private)
        videoInput;
        videoSource;
    end
    properties (Dependent = true)
        integrationTime; % seconds
        gain;
        regionOfInterest;
    end
    
    methods
        function cam=PikeCam(deviceID)
            cam=cam@Cam();
            if (nargin<1 || isempty(deviceID))
                deviceID=-1;
            end
            
            %Check if there is already a GigE cam configured
            vids=imaqfind('Tag','PikeVideoinput');
            vid=[];
            vidIdx=0;
            while (vidIdx<length(vids) && isempty(vid))
                vidIdx=vidIdx+1;
                vid=vids{vidIdx};
                if (deviceID>0 && vid.DeviceID~=deviceID)
                    vid=[];
                end
            end
            if (isempty(vid))
                hwInfo=imaqhwinfo('avtmatlabadaptor64_r2009b');
                if (deviceID<0 && ~isempty(hwInfo.DeviceIDs))
                    deviceID=hwInfo.DeviceIDs{1};
                end
                deviceInfoIdx=find(cell2mat(hwInfo.DeviceIDs)==deviceID);
                if (~isempty(deviceInfoIdx))
                    %formats=hwInfo.DeviceInfo(deviceInfoIdx(1)).SupportedFormats;
                    vid=videoinput('avtmatlabadaptor64_r2009b',deviceID,'F0M6_Mono16_640x480');
                    vid.Tag='PikeVideoinput';
                else
                    error('No Pike camera found. Please check if the driver is correctly installed.');
                end
            end
            %Retrieve the number of bits per pixel from the format string
            tokens=regexp(get(vid,'VideoFormat'),'Mono([\d]+)_','tokens');
            cam.bitsPerPixel=str2double(tokens{1});
            
            src=getselectedsource(vid);
            
            src.Shutter=500; %[0:4095]
            src.ExtendedShutter=50000; % [0:((2^26)-1)]
            src.Gain=0;
            src.Brightness=16;
            src.Timebase=4;
                        
            cam.videoInput=vid;
            cam.videoSource=src;
            
            ensureVideoStopped(cam);
            
            cam.videoInput.ROIPosition=[0 0 cam.videoInput.VideoResolution];
            exposureTime=5e-3; %In seconds
            cam.videoSource.ExtendedShutter=1e6*exposureTime;
            
            set(cam.videoInput,'TriggerRepeat',Inf);
            triggerconfig(cam.videoInput,'manual');
            ensureVideoStarted(cam);
            
            cam.maxSize([2 1])=cam.videoInput.VideoResolution;
            
            cam.regionOfInterest=[0 0 cam.maxSize];
            
            cam.numberOfFramesToAverage=1;
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            props=propinfo(cam.videoSource,'ExtendedShutter');
            if (~isempty(props.ConstraintValue))
                integrationTimeLimits=1e-6*(props.ConstraintValue-[0 1]);
            else
                integrationTimeLimits=[1e-9 10];
            end
            newIntegrationTime=min(integrationTimeLimits(2),max(integrationTimeLimits(1),newIntegrationTime));
            
            ensureVideoStopped(cam);
            cam.videoSource.ExtendedShutter=1e6*newIntegrationTime;
            ensureVideoStarted(cam);
        end
        function integrationTime=get.integrationTime(cam)
            integrationTime=cam.videoSource.ExtendedShutter*1e-6;
        end
        function cam=set.gain(cam,newGain)
            props=propinfo(cam.videoSource,'Gain');
            gainLimits=props.ConstraintValue;
            if (isempty(gainLimits))
                gainLimits=[0 630];
            end
            newGain=min(gainLimits(2),max(gainLimits(1),newGain));
            
            ensureVideoStopped(cam);
            cam.videoSource.Gain=newGain;
            ensureVideoStarted(cam);
        end
        function gain=get.gain(cam)
            gain=cam.videoSource.Gain;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=[];
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            ensureVideoStopped(cam);
            %Make sure that the size is at least 128x128
            sizeIncrease=max(0,cam.minRegionOfInterestSize-newROI(3:4));
            newROI(1:2)=newROI(1:2)+floor(-sizeIncrease/2);
            newROI(3:4)=newROI(3:4)+sizeIncrease;
            %Round to even numbers when using a Bayer filter
            newROI=2*floor(newROI./2);
            %Clip the region of interest to the active area of the CCD
            newROI(1:2)=min(max(0,newROI(1:2)),cam.maxSize-1);
            newROI(3:4)=min(newROI(1:2)+max(2,newROI(3:4)),cam.maxSize)-newROI(1:2);
            cam.videoInput.ROIPosition=newROI([2 1 4 3]);
            ensureVideoStarted(cam);
        end
        function roi=get.regionOfInterest(cam)
            roi([2 1 4 3])=cam.videoInput.ROIPosition;
        end
        function delete(cam)
            delete@Cam(cam);
            
            ensureVideoStopped(cam);
            delete(cam.videoInput);
            clear cam.videoInput;
        end
    end
    methods(Access = protected)
        function img=acquireDirect(cam,nbFrames)
            nbColorChannels=1;
            imgSize=cam.regionOfInterest(3:4);
            img=zeros([imgSize nbColorChannels nbFrames]);
            
            % The following is not thread-safe, so better store this now
            numberOfFramesToAverage=cam.numberOfFramesToAverage;
            for (frameIdx=1:nbFrames)
                for (triggerIdx=1:cam.numberOfFramesToAverage)
                    snapshot=[];
                    while (isempty(snapshot) || size(snapshot,1)~=imgSize(1) || size(snapshot,2)~=imgSize(2))
                        ensureVideoStarted(cam);
                        try
                            snapshot=single(swapbytes(getsnapshot(cam.videoInput)));
                        catch Exc
                            logMessage('Something went wrong, reseting the videoinput object...');
                            stop(cam.videoInput);
                            start(cam.videoInput);
                        end
                    end
                    %Normalize the maximum graylevel to 1
                    snapshot=snapshot./(2^cam.bitsPerPixel-4);

                    %Debug info
                    if (max(snapshot(:))==1)
                        logMessage('%u pixels are saturated in %u frames!',[sum(snapshot(:)==max(snapshot(:))) framesPerTrigger]);
                    end

                    img(:,:,:,frameIdx)=img(:,:,:,frameIdx)+sum(snapshot,4);
                end
            end
            
            img=img./numberOfFramesToAverage;
        end
    end
    methods (Access = private)
        function ensureVideoStarted(cam)
            if (~strcmpi(cam.videoInput.Running,'on'))
                start(cam.videoInput);
            end
        end
        function ensureVideoStopped(cam)
            if (strcmpi(cam.videoInput.Running,'on'))
                stop(cam.videoInput);
            end
        end
    end
    
end

