%
classdef AndorCam < Cam
    %Orca Flash 4 camera class
    %  
    properties (SetAccess = private)
        maxSize;
        bitsPerPixel;
    end
    properties (Access = private)
       hndl;
       rc;
       stride;
       imagesize; %length of buffer
       running = false
       softwareRegionOfInterest;
       minRegionOfInterestSize=[1 1]*128;
       hardwareRegionOfInterest;
       hardWareRegionOfInterestStepSize=[1 1]*4;
    end
    properties (Dependent = true)
        integrationTime; % seconds
        gain;
        regionOfInterest;
    end
 
    
    
    methods
        function cam=AndorCam(deviceID)
            cam=cam@Cam()
            if (nargin<1 || isempty(deviceID))
                deviceID=0;
            end
            
            %Initialise camera
            [rc] = AT_InitialiseLibrary();
            AT_CheckError(rc);
            [rc,cam.hndl] = AT_Open(deviceID); % default 0
            AT_CheckError(rc);
            
            ensureVideoStopped(cam);
            
            [rc] = AT_SetEnumString(cam.hndl,'CycleMode','Continuous');
            AT_CheckWarning(rc);
            [rc] = AT_SetEnumString(cam.hndl,'TriggerMode','Software');
            AT_CheckWarning(rc);
            disp('Camera initialized');
        
            %set the number of bits per pixel to 16
            cam.bitsPerPixel=16;
            [rc] = AT_SetEnumString(cam.hndl,'SimplePreAmpGainControl','16-bit (low noise & high well capacity)');
            AT_CheckWarning(rc);
            [rc] = AT_SetEnumString(cam.hndl,'PixelEncoding','Mono16');
            AT_CheckWarning(rc);            
            
            ensureVideoStarted(cam);
            
            %retrival parameters from camera

            cam.maxSize = getcurrentSize(cam);
            cam.regionOfInterest=[0 0 cam.maxSize];
            [rc,cam.imagesize] = AT_GetInt(cam.hndl,'ImageSizeBytes');
            AT_CheckWarning(rc);
            [rc,cam.stride] = AT_GetInt(cam.hndl,'AOIStride'); 
            AT_CheckWarning(rc);
            %
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            ensureVideoStopped(cam);
            
            newIntegrationTime = 0.02;
            [rc] = AT_SetFloat(cam.hndl,'ExposureTime',newIntegrationTime);
            AT_CheckWarning(rc);
            ensureVideoStarted(cam);
        end    
        function integrationTime=get.integrationTime(cam)
            [rc, integrationTime] = AT_GetFloat(cam.hndl,'ExposureTime');
            AT_CheckWarning(rc);
        end
        function cam=set.gain(cam,newGain)
            %TODO
        end
        function gain=get.gain(cam)
            gain=1;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=[];
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            %Store this value for later software clipping or expanding
            cam.softwareRegionOfInterest=newROI;
            
            %Clip the region of interest to the active area of the sCMOS
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
            
            %Set up AOI 
%            ensureVideoStopped(cam);            
%             [rc] = AT_SetInt(cam.hndl, 'AOIHeight',newROI(3));
%             AT_CheckWarning(rc);
%             [rc] = AT_SetInt(cam.hndl, 'AOIWidth',newROI(4));
%             AT_CheckWarning(rc);
%             [rc] = AT_SetInt(cam.hndl, 'AOITop',newROI(1));
%             AT_CheckWarning(rc);
%             [rc] = AT_SetInt(cam.hndl, 'AOILeft',newROI(2));
%             AT_CheckWarning(rc);
%            ensureVideoStarted(cam);
            
            
        end
        
        function roi=get.regionOfInterest(cam)
            roi=cam.softwareRegionOfInterest;
        end
        
        function delete(cam)
            ensureVideoStopped(cam);
            delete@Cam(cam);    
            [rc] = AT_Flush(cam.hndl);
            AT_CheckWarning(rc);
            [rc] = AT_Close(cam.hndl);
            AT_CheckWarning(rc);
            [rc] = AT_FinaliseLibrary();
            AT_CheckWarning(rc);
            disp('Camera shutdown');
        end
        function size=getcurrentSize(cam)
            [rc,height] = AT_GetInt(cam.hndl,'AOIHeight');
            AT_CheckWarning(rc);
            [rc,width] = AT_GetInt(cam.hndl,'AOIWidth');  
            AT_CheckWarning(rc);
            size=[height, width];
        end
    end
    
    
    methods(Access = protected)
        function img=acquireDirect(cam,nbFrames)
            
            if (nargin<2 || isempty(nbFrames))
                nbFrames = 1;
            end           
            maxPixelValue=(2^cam.bitsPerPixel-1);
            size = getcurrentSize(cam);
            ensureVideoStarted(cam);
            %img = zeros(width,height);
            for X = 1:nbFrames               
                [rc] = AT_QueueBuffer(cam.hndl,cam.imagesize);
                AT_CheckWarning(rc);
                [rc] = AT_Command(cam.hndl,'SoftwareTrigger');
                AT_CheckWarning(rc);
                [rc,buf] = AT_WaitBuffer(cam.hndl,1000);
                AT_CheckWarning(rc);
                [rc,img] = AT_ConvertMono16ToMatrix(buf,size(1),size(2),cam.stride);
            end
            img=double(img)./maxPixelValue;
        end
    end
    
    
    
    methods (Access = private)

        function ensureVideoStarted(cam)
            if (cam.running == false)
                [rc] = AT_Command(cam.hndl,'AcquisitionStart');
                AT_CheckWarning(rc);
                cam.running = true;
            end
        end
        function ensureVideoStopped(cam)
            if (cam.running ~= false)
                [rc] = AT_Command(cam.hndl,'AcquisitionStop');
                AT_CheckWarning(rc);
                cam.running = false;
            end
        end
            
    end
end