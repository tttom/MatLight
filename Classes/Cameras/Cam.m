classdef Cam < handle
    % Camera super class
    %
    properties (Abstract)
        integrationTime;
        gain;
        regionOfInterest;
    end
    properties (Abstract, SetAccess = private)
        maxSize;
        bitsPerPixel;
    end
    properties
        background=[];
        acquisitionFunctor=[];
        numberOfFramesToAverage=1;
    end
    
    methods
        function cam=Cam()
        end
        function cam=acquireBackground(cam,nbFramesToAverage)
            if (nargin<2)
                nbFramesToAverage=cam.numberOfFramesToAverage;
                % Temporarily change the number of frames that needs to be averaged
                origNbFramesToAverage=cam.numberOfFramesToAverage;
                cam.numberOfFramesToAverage=nbFramesToAverage;
            end
            cam.background = cam.acquireDirect(1);
            if (nargin<2)
                cam.numberOfFramesToAverage=origNbFramesToAverage;
            end
        end
        function img=acquire(cam,nbFrames)
            if (nargin<2)
                nbFrames=1;
            end
            img=cam.acquireDirect(nbFrames);
            if (~isempty(cam.background))
                img=img-repmat(cam.background,[1 1 1 nbFrames]);
            end
            if (~isempty(cam.acquisitionFunctor))
                if (ischar(cam.acquisitionFunctor))
                    cam.acquisitionFunctor=str2func(cam.acquisitionFunctor);
                end
                switch (nargin(cam.acquisitionFunctor))
                    case 0
                        cam.acquisitionFunctor();
                    otherwise
                        cam.acquisitionFunctor(img);
                end
            end
        end
    end
    methods (Abstract=true, Access = protected)
        img=acquireDirect(cam,nbFrames);
    end 
end