%
% Mikrotron camera class
% 
classdef MikrotronCam < Cam
    properties (SetAccess = private)
        maxSize;
        minRegionOfInterestSize=[1 1]*2;
        bitsPerPixel;
    end
    properties (Access = private)
        videoInput;
        videoSource;
        canSetIntegrationTime;
        canSetGain;
    end
    properties (Dependent = true)
        integrationTime; % seconds
        gain;
        regionOfInterest;
    end
    
    methods
        function cam=MikrotronCam(deviceID)
            cam=cam@Cam();
            if (nargin<1 || isempty(deviceID))
                deviceID=-1;
            end
            
            %Check if there is already a GigE cam configured
            vids=imaqfind('Tag','MikrotronVideoinput');
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
                hwInfo=imaqhwinfo('ni');
                if (deviceID<0 && ~isempty(hwInfo.DeviceIDs))
                    deviceID=hwInfo.DeviceIDs{1};
                end
                deviceInfoIdx=find(cell2mat(hwInfo.DeviceIDs)==deviceID);
                if (~isempty(deviceInfoIdx))
                    formats=hwInfo.DeviceInfo(deviceInfoIdx(1)).SupportedFormats;
                    if (any(strcmp(formats,'ni')))
                        vid=videoinput('ni',deviceID,'img0');
                        vid.ReturnedColorspace = 'grayscale';
                        %vid.BayerSensorAlignment = 'grbg';                   
                    else
                        devInfo=imaqhwinfo('ni',deviceID);
                        vid=videoinput('ni',deviceID,devInfo.DefaultFormat);
                        set(vid,'ReturnedColorSpace','grayscale');
                    end
                    vid.Tag='MikrotronVideoinput';
                else
                     error('No Mikrotron camera found.');
                end
            end
            %Retrieve the number of bits per pixel from the format string
            cam.bitsPerPixel=12;
            
            src=getselectedsource(vid);
            
            cam.videoInput=vid;
            cam.videoSource=src;
            
            ensureVideoStopped(cam);
            
            cam.videoInput.ROIPosition=[0 0 cam.videoInput.VideoResolution];
                          
            set(cam.videoInput,'TriggerRepeat',Inf);
            triggerconfig(cam.videoInput,'manual');
            ensureVideoStarted(cam);
            
            cam.maxSize([2 1])=cam.videoInput.VideoResolution;
            
            cam.regionOfInterest=[0 0 cam.maxSize];
            cam.canSetIntegrationTime = 0;
            cam.canSetGain = 0;
            cam.numberOfFramesToAverage=1;
        end
        
        function cam=set.integrationTime(cam,newIntegrationTime)
            logMessage('Could not set exposure, Mikrotron driver cannot supput it.');
            cam.integrationTime=-1;
        end
        function integrationTime=get.integrationTime(cam)
            logMessage('Could not set exposure, Mikrotron driver cannot supput it.');
            integrationTime=-1;
        end
        function cam=set.gain(cam,newGain)
            if (cam.canSetGain)
                ensureVideoStopped(cam);
                registerValue=log10(newGain)*1023/1.8; % 10 bits value linear with gain in dB with maximum of 18dB it seems
                % 36dB accorduing to http://www.theimagingsourceforums.com/archive/index.php/t-324235.html
                cam.videoSource.Gain=max(260,min(1023,round(registerValue))); %260-1023;
                ensureVideoStarted(cam);
            end
        end
        function gain=get.gain(cam)
            if (cam.canSetGain)
                gain=10^(1.8*cam.videoSource.Gain/1023);
            else
                gain=0;
            end
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
            maxPixelValue=(2^cam.bitsPerPixel-1-1.0*cam.maxValue254);
            imgSize=cam.hardwareRegionOfInterest(3:4);
            framesPerTrigger=nbFrames*cam.numberOfFramesToAverage;
             
            img=[];
            while (isempty(img) || any([size(img,1) size(img,2)]~=imgSize))
                ensureVideoStarted(cam);
                trigger(cam.videoInput);
                try
                    img=getdata(cam.videoInput,framesPerTrigger);
                    %Debug info
                    if (max(img(:))>=maxPixelValue)
                        msg=sprintf('%u pixels are saturated',sum(img(:)==max(img(:))));
                        if (framesPerTrigger~=1)
                            msg=strcat(msg,sprintf(' in %u frames',framesPerTrigger));
                        end
                        logMessage(msg);
                    end
                    % Average cam.numberOfFramesToAverage consecutive frames
                    img=mean(single(reshape(img,[],cam.numberOfFramesToAverage,nbFrames)),2);
                    img=reshape(img,[imgSize nbColorChannels nbFrames]);
                catch Exc
                    logMessage('Something went wrong, reseting the videoinput object...');
                    stop(cam.videoInput);
                    start(cam.videoInput);
                end
            end
            %Normalize the maximum graylevel to 1
            img=img./maxPixelValue;
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


