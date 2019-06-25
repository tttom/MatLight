classdef DirectShowCam < Cam
    %DirectShow camera class
    %Uses vcapg2.dll to access DirectShow
    %  
    properties
        integrationTime = 20e-3; % seconds
        gain;
        regionOfInterest;
    end
    properties (SetAccess = private)
        maxSize;
        validRegionOfInterest;
        bitsPerPixel;
    end
    properties (Access = private)
        cardNumber;
    end
    
    methods
        function cam=DirectShowCam()
            cam=cam@Cam();
            if (exist('vcapg2'))
                cam.cardNumber=vcapg2();
                testImg=vcapg2(1);
                cam.maxSize=[size(testImg,1) size(testImg,2)];
                cam.regionOfInterest=[0 0 cam.maxSize];
                cam.numberOfFramesToAverage=1;
                cam.bitsPerPixel=8;
            else
                logMessage('Please install vcapg2 from http://www.mathworks.com/matlabcentral/fileexchange/2939');
            end
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            logMessage('vcapg2 doesn''t permit changing the integration time.');
            cam.integrationTime=newIntegrationTime;
        end
        function integrationTime=get.integrationTime(cam)
            logMessage('vcapg2 doesn''t permit querying the integration time.');
            integrationTime=cam.integrationTime;
        end
        function cam=set.gain(cam,newGain)
            logMessage('vcapg2 doesn''t permit changing the gain.');
            cam.gain=newGain;
        end
        function gain=get.gain(cam)
            logMessage('vcapg2 doesn''t permit querying the gain.');
            gain=cam.gain;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=[];
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            newROI(1:2)=min(cam.maxSize,max(0,newROI(1:2)));
            newROI(3:4)=min(cam.maxSize,max(1,newROI(3:4)+newROI(1:2)))-newROI(1:2);
            cam.validRegionOfInterest=newROI;
        end
        function roi=get.regionOfInterest(cam)
            roi=cam.validRegionOfInterest;
        end
        function delete(cam)
            delete@Cam(cam);
        end
    end
    methods (Access = protected)
        function img=acquireDirect(cam,nbFrames)
            maxPixelValue=2^cam.bitsPerPixel-1;
            nbColorChannels=1;
            imgSize=cam.regionOfInterest(3:4);
            img=zeros([imgSize nbColorChannels nbFrames]);
             
            % The following is not thread-safe, so better store this now
            numberOfFramesToAverage=cam.numberOfFramesToAverage;
            regionOfInterest=cam.regionOfInterest;
            for (frameIdx=1:nbFrames),
                for (triggerIdx=1:cam.numberOfFramesToAverage)
                    snapshot=[];
                    nbTrials=5;
                    while (isempty(snapshot) || (nbTrials>0 && all(snapshot(:)==0)))
                        try
                            snapshot=vcapg2(1);
                        catch Exc
                            logMessage('Something went wrong, trying again...');
                        end
                        nbTrials=nbTrials-1;
                    end
                    if (nbTrials==0)
                        logMEssage('Received all zero image from camera.');
                    end
                    %Select the region of interest
                    snapshot=single(snapshot(regionOfInterest(1)+[1:regionOfInterest(3)],regionOfInterest(2)+[1:regionOfInterest(4)],:));

                    %Debug info
                    if (max(snapshot(:))>=maxPixelValue)
                        logMessage('%u pixels are saturated!',sum(snapshot(:)==max(snapshot(:)))./3);
                    end
                    
                    snapshot=0.2126*snapshot(:,:,1)+0.7152*snapshot(:,:,2)+0.0722*snapshot(:,:,3); %Objective luminance weighting

                    img(:,:,:,frameIdx)=img(:,:,:,frameIdx)+snapshot;
                end
            end
            img=img./numberOfFramesToAverage;
            img=img./maxPixelValue;
        end
    end
    
end