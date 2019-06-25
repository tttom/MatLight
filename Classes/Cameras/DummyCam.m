classdef DummyCam < Cam
    %Dummy camera class
    %
    properties
        integrationTime = 20e-3;
        gain = 1;
        regionOfInterest;
        acquireDirectFunctor=[];
    end
    properties (SetAccess = private)
        maxSize = [256 256];
        bitsPerPixel = 8;
    end
    properties (Access = private)
        testObject;
    end
    
    methods
        function cam=DummyCam(imgSize)
            cam=cam@Cam();
            
            if (nargin<1 || isempty(imgSize))
                imgSize=[1 1]*128;
            end
            
            cam.maxSize=imgSize;
            
            cam.regionOfInterest=[0 0 cam.maxSize];
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            cam.integrationTime=max(1e-9,min(60,newIntegrationTime));
        end
        function integrationTime=get.integrationTime(cam)
            integrationTime=cam.integrationTime;
        end
        function cam=set.gain(cam,newGain)
            cam.gain=max(10^-3,min(10^6,newGain));
        end
        function gain=get.gain(cam)
            gain=cam.gain;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=[];
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            cam.regionOfInterest=newROI;
            
            %Adjust the test object
            cam.testObject=zeros(cam.maxSize);
            cam.testObject(floor(cam.maxSize(1)/2)+1,:)=100;
            cam.testObject(:,floor(cam.maxSize(2)/2)+1)=100;
            cam.testObject=cam.testObject(cam.regionOfInterest(1)+[1:cam.regionOfInterest(3)],cam.regionOfInterest(2)+[1:cam.regionOfInterest(4)]);
        end
        function cam=acquireBackground(cam,nbFrames)
            if (nargin<2)
                nbFrames=cam.numberOfFramesToAverage;
            end
            
            %Remove the test object before taking the dark image
            testObject=cam.testObject;
            cam.testObject=0*testObject;
            cam=acquireBackground@Cam(cam,nbFrames);
            cam.testObject=testObject;
        end
        function img=acquire(cam,nbFrames)
            if (nargin<2)
                nbFrames=1;
            end
            img=acquire@Cam(cam,nbFrames);
        end
        function delete(cam)
            delete@Cam(cam);
        end
    end
    methods (Access = protected)
        function img=acquireDirect(cam,nbFrames)
            if (isempty(cam.acquireDirectFunctor))
                img=cam.acquireDirectDefaultFunctor(nbFrames,cam.numberOfFramesToAverage);
            else
                switch nargin(cam.acquireDirectFunctor)
                    case 0
                        img=cam.acquireDirectFunctor();
                    case 1
                        img=cam.acquireDirectFunctor(cam);
                    case 2
                        img=cam.acquireDirectFunctor(cam,nbFrames);
                    otherwise
                        img=cam.acquireDirectFunctor(cam,nbFrames,cam.numberOfFramesToAverage);
                end
            end
        end
    end
    methods (Access = private)
        function img=acquireDirectDefaultFunctor(cam,nbFrames,numberOfFramesToAverage)
            integrationTimeUnits=nbFrames*cam.integrationTime/20e-3;
            
            frm=repmat(cam.testObject,[1 1 1 nbFrames]);
            frm=frm+5; %Dark noise
            frm=frm*cam.gain;
            frm=frm*integrationTimeUnits+sqrt(frm*integrationTimeUnits).*randn(size(frm))./sqrt(numberOfFramesToAverage);
            
            img=floor(frm);
            
            %Normalize the maximum graylevel to 1
            img=img./(2^cam.bitsPerPixel-1);
        end
    end
end

