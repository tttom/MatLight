%
classdef BaslerGigECam < Cam
    %Basler GigE camera class
    %  
    properties (SetAccess = private)
        maxSize;
        minRegionOfInterestSize=[1 1]*128;
        bitsPerPixel;
    end
    properties (Access = private)
        videoInput;
        videoSource;
        softwareRegionOfInterest;
        hardwareRegionOfInterest;
        maxValue254=false;
    end
    properties (Dependent = true)
        integrationTime; % seconds
        gain;
        regionOfInterest;
    end
    
    methods
        function cam=BaslerGigECam(deviceID)
            cam=cam@Cam();
            if (nargin<1 || isempty(deviceID))
                deviceID=-1;
            end
            
            %Check if there is already a GigE cam configured
            vids=imaqfind('Tag','GigEVideoinput');
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
                hwInfo=imaqhwinfo('gige');
                if (deviceID<0 && ~isempty(hwInfo.DeviceIDs))
                    deviceID=hwInfo.DeviceIDs{1};
                end
                deviceInfoIdx=find(cell2mat(hwInfo.DeviceIDs)==deviceID);
                if (~isempty(deviceInfoIdx))
                    formats=hwInfo.DeviceInfo(deviceInfoIdx(1)).SupportedFormats;
                    if (any(strcmp(formats,'Mono12Packed')))
                        vid=videoinput('gige',deviceID,'Mono12Packed');
                    else
                        vid=videoinput('gige',deviceID,'Mono8');
                    end
                    if (strcmpi(hwInfo.DeviceInfo(deviceInfoIdx(1)).DeviceName(end),'c'))
                        cam.maxValue254=true;
                    else
                        cam.maxValue254=false;
                    end
                    vid.Tag='GigEVideoinput';
                else
                    error('No GigE camera found. Please check if: it is plugged in power and ethernet, the (Pylon) filter driver is unchecked in network properties, it is installed in Matlab using installgenicam.');
                end
            end
            %Retrieve the number of bits per pixel from the format string
            cam.bitsPerPixel=str2double(regexprep(get(vid,'VideoFormat'),'[^\d]',''));
            
            src=getselectedsource(vid);
            
            src.TestImageSelector='Off';
            
            cam.videoInput=vid;
            cam.videoSource=src;
            
            ensureVideoStopped(cam);
            
            cam.videoInput.ROIPosition=[0 0 cam.videoInput.VideoResolution];
            exposureTime=5e-3; %In seconds
            cam.videoSource.ExposureTimeAbs=1e6*exposureTime;
            cam.videoSource.AllGainRaw=0; %130;
            if any(strcmp(fieldnames(cam.videoSource),'ExposureAuto'))
                cam.videoSource.ExposureAuto = 'Off';
            end
            if any(strcmp(fieldnames(cam.videoSource),'GainAuto'))
                cam.videoSource.GainAuto = 'Off';
            end
            if any(strcmp(fieldnames(cam.videoSource),'PacketSize'))
                if ischar(cam.videoSource.PacketSize)
                    cam.videoSource.PacketSize ='9014';
                else
                    cam.videoSource.PacketSize=8192;
                end
            end
            set(cam.videoInput,'TriggerRepeat',Inf);
            triggerconfig(cam.videoInput,'manual');
            ensureVideoStarted(cam);
            
            cam.maxSize([2 1])=cam.videoInput.VideoResolution;
            
            cam.regionOfInterest=[0 0 cam.maxSize];
            
            cam.numberOfFramesToAverage=1;
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            props=propinfo(cam.videoSource,'ExposureTimeAbs');
            if (~isempty(props.ConstraintValue))
                integrationTimeLimits=1e-6*(props.ConstraintValue-[0 1]);
            else
                integrationTimeLimits=[1e-9 10];
            end
            newIntegrationTime=min(integrationTimeLimits(2),max(integrationTimeLimits(1),newIntegrationTime));
            
            ensureVideoStopped(cam);
            cam.videoSource.ExposureTimeAbs=1e6*newIntegrationTime;
            ensureVideoStarted(cam);
        end
        function integrationTime=get.integrationTime(cam)
            integrationTime=cam.videoSource.ExposureTimeAbs*1e-6;
        end
        function cam=set.gain(cam,newGain)
            props=propinfo(cam.videoSource,'AllGainRaw');
            gainLimits=props.ConstraintValue;
            if (isempty(gainLimits))
                gainLimits=[0 500];
            end
            newGain=min(gainLimits(2),max(gainLimits(1),newGain));
            
            ensureVideoStopped(cam);
            cam.videoSource.AllGainRaw=newGain;
            ensureVideoStarted(cam);
        end
        function gain=get.gain(cam)
            gain=cam.videoSource.AllGainRaw;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=[];
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            %Store this value for later software clipping or expanding
            cam.softwareRegionOfInterest=newROI;
            
            %Round to even numbers for when using a Bayer filter
            newROI=2*[floor(newROI(1:2)./2) ceil(floor(newROI(3:4)./2))];
            
            %Clip the region of interest to the active area of the CCD
            newROI(1:2)=min(max(0,newROI(1:2)),cam.maxSize-1);
            newROI(3:4)=min(newROI(1:2)+max(2,newROI(3:4)),cam.maxSize)-newROI(1:2);
            
            %Make sure that the size is at least minRegionOfInterestSize
            if (any(newROI(3:4)<cam.minRegionOfInterestSize))
                sizeIncrease=max(0,cam.minRegionOfInterestSize-newROI(3:4));
                newROI(3:4)=newROI(3:4)+sizeIncrease;
                %Make sure the ROI fits within the detector range
                newROI(1:2)=min(cam.maxSize-newROI(3:4),max(0,newROI(1:2)+floor(-sizeIncrease/2)));
            end
            cam.hardwareRegionOfInterest=newROI;
                        
            % Update the camera driver ROI
            ensureVideoStopped(cam);
            cam.videoInput.ROIPosition=newROI([2 1 4 3]);
            ensureVideoStarted(cam);
        end
        function roi=get.regionOfInterest(cam)
            roi=cam.softwareRegionOfInterest;
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
             
            % The following is not thread-safe, so better store this now
            numberOfFramesToAverage=cam.numberOfFramesToAverage;
            framesPerTrigger=nbFrames*numberOfFramesToAverage;
            if (get(cam.videoInput,'FramesPerTrigger')~=framesPerTrigger)
                ensureVideoStopped(cam);
                set(cam.videoInput,'FramesPerTrigger',framesPerTrigger);
            end
                
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
                    img=mean(single(reshape(img,[],numberOfFramesToAverage,nbFrames)),2);
                    img=reshape(img,[imgSize nbColorChannels nbFrames]);
                catch Exc
                    logMessage('Something went wrong, reseting the videoinput object...');
                    stop(cam.videoInput);
                    start(cam.videoInput);
                end
            end
            %Normalize the maximum graylevel to 1
            img=img./maxPixelValue;
            
            %Check if we need to crop or expand the image in software
            if (any(cam.softwareRegionOfInterest~=cam.hardwareRegionOfInterest))
                rangeBeginDiffs=cam.softwareRegionOfInterest(1:2)-cam.hardwareRegionOfInterest(1:2);
                rangeEndDiffs=cam.softwareRegionOfInterest(1:2)+cam.softwareRegionOfInterest(3:4)-cam.hardwareRegionOfInterest(1:2)-cam.hardwareRegionOfInterest(3:4);
                % Crop
                img=img(max(1,1+max(0,rangeBeginDiffs(1))):end-max(0,-rangeEndDiffs(1)),max(1,1+max(0,rangeBeginDiffs(2))):end-max(0,-rangeEndDiffs(2)),:,:);
                if (any([size(img,1) size(img,2)]<cam.softwareRegionOfInterest(3:4)))
                    % Zero pad
                    resizedImg=zeros([cam.softwareRegionOfInterest(3:4) size(img,3) size(img,4)]);
                    resizedImg(max(0,-rangeBeginDiffs(1))+[1:size(img,1)],max(0,-rangeBeginDiffs(2))+[1:size(img,2)],:,:)=img;
                    img=resizedImg; clear resizedImg;
                end
            end
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

