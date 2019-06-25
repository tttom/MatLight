% Note, use the following to close all full-screens:
%     CogentGraphicsScreenManager.instance().delete();
classdef CogentGraphicsScreenManager < ScreenManager
    properties (SetAccess=private, Dependent=true)
        screens;
    end
    properties (Access=private)
        currentScreen;
    end
    methods(Access=private)
        function sm=CogentGraphicsScreenManager()
            sm.currentScreen=[];
        end
    end
    methods (Static)
        function obj = instance()
            persistent uniqueInstance;
            if isempty(uniqueInstance)
                obj = CogentGraphicsScreenManager();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    methods
        function scrns=get.screens(sm)
            scrns=struct([]);
            ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
            screenDevices = ge.getScreenDevices();
            defaultScreen = ge.getDefaultScreenDevice();
            for (displayIdx=1:length(screenDevices))
                displaySize(1) = screenDevices(displayIdx).getDisplayMode().getHeight();
                displaySize(2) = screenDevices(displayIdx).getDisplayMode().getWidth();
                scrns(displayIdx).size=displaySize;
                scrns(displayIdx).frameRate=screenDevices(displayIdx).getDisplayMode().getRefreshRate();
                scrns(displayIdx).open=sm.isOpen(displayIdx);
                scrns(displayIdx).id=displayIdx;
                scrns(displayIdx).main=screenDevices(displayIdx)==defaultScreen;
                scrns(displayIdx).description=sprintf('%u (%ux%u)',[displayIdx displaySize([2 1])]);
                if (scrns(displayIdx).main)
                    scrns(displayIdx).description=[scrns(displayIdx).description,' main'];
                end
            end
        end
        function sm=display(sm,id,img,position)
            if (nargin==1)
                % Special case Matlab built-in: just display a text string
                scrnsStr='';
                for id=1:length(sm.screens);
                    scr=sm.screens(id);
                    if ~isempty(scr)
                        scrnsStr=strcat(scrnsStr,', ',scr.description);
                    end
                end
                scrnsStr=scrnsStr(3:end);
                description=strcat('CogentGraphicsScreenManager: ',scrnsStr);
                disp(description);
                methods(sm);
                properties(sm);
                return;
            end
            if (~sm.isOpen(id))
                sm=sm.open(id);
            end
            if (nargin<4 || isempty(position))
                position=[0 0 sm.currentScreen.size];
            end
                        
            %Crop image
            maxRoiSize=sm.currentScreen.size-position(1:2);
            if (size(img,1)>maxRoiSize(1))
                img=img(1:maxRoiSize(1),:,:);
            end
            if (size(img,2)>maxRoiSize(2))
                img=img(:,1:maxRoiSize(2),:);
            end
            %Convert to true color if in grayscale.
            if (size(img,3)>1),
                trueColorImageMatrix=img;
            else
                trueColorImageMatrix=repmat(img,[1 1 3]);%Slow and memory inefficient
            end
            
            imageHeight=size(img,1);
            imageWidth=size(img,2);
            cgloadarray(1,imageWidth,imageHeight,reshape(permute(trueColorImageMatrix,[2 1 3]),imageWidth*imageHeight,[]),imageWidth,imageHeight);
            cgdrawsprite(1,position(2)-floor((sm.currentScreen.size(2)-imageWidth)./2),floor((sm.currentScreen.size(1)-imageHeight)./2)-position(1)); %Careful! Odd convention for X direction!
            cgflip(0,0,0); %Also blank out next frame
        end
    end
    methods (Access=protected)
        function sm=openSingle(sm,id)
            if (~sm.isOpen(id))
                availableScreens=sm.screens;
                if (id<=length(availableScreens))
                    %Close any already open screen first because CG can only
                    %handle one screen at a time
                    if (~isempty(sm.currentScreen) && sm.isOpen(sm.currentScreen.id))
                        sm.close(sm.currentScreen.id);
                    end
                    % Open the new screen
                    screenSize=sm.screens(id).size;
                    frameRate=sm.screens(id).frameRate; %0=fastest, or 56-60; % Hz

                    cgloadlib();                
                    cgopen(screenSize(2),screenSize(1),0,frameRate,id);

                    sm.currentScreen.id=id;
                    sm.currentScreen.size=screenSize;
                else
                    error('Could not detect screen %d, unable to open.',id);
                end
            else
                logMessage('Full screen %d is already opened',id);
            end
        end
        function sm=closeSingle(sm,id)
            if (sm.isOpen(id))
                try
                    cgshut();
                    sm.currentScreen=[];
                catch Exc
                    logMessage('Could not close full screen window.');
                    rethrow(Exc);
                end
            else
                logMessage('Screen not open.');
            end
            
        end
    end
    methods (Access=private)
        function TF=isOpen(sm,id)
            TF=~isempty(sm.currentScreen) && sm.currentScreen.id==id;
        end
    end
end