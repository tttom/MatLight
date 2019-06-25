% Note, use the following to close all full-screens:
%     JavaScreenManager.instance().delete();
classdef JavaScreenManager < ScreenManager
    properties (SetAccess=private, Dependent=true)
        screens;
    end
    properties (Access=private)
        fullScreenWindows;
    end
    methods(Access=private)
        function sm=JavaScreenManager()
            sm.fullScreenWindows=struct([]);
        end
    end
    methods (Static)
        function obj = instance()
            persistent uniqueInstance;
            if isempty(uniqueInstance)
                obj = JavaScreenManager();
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
                description=strcat('JavaScreenManager: ',scrnsStr);
                disp(description);
                methods(sm);
                properties(sm);
                return;
            end
            if (~sm.isOpen(id))
                sm=sm.open(id);
            end
            if (nargin<4 || isempty(position))
                position=[0 0 sm.fullScreenWindows(id).size];
            end
            fullScreenBufferedImage=sm.fullScreenWindows(id).bufferedImage;
                        
            %Crop image
            maxRoiSize=sm.fullScreenWindows(id).size-position(1:2);
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
            
            rgbaImg=trueColorImageMatrix;
            rgbaImg=min(255,floor(256*rgbaImg));
%             rgbaImg=255*256^3+reshape(rgbaImg,[],3)*(256.^[2 1 0].'); % alpha channel
            rgbaImg=reshape(rgbaImg,[],3)*(256.^[2 1 0].');
            rgbaImg=reshape(rgbaImg,size(img,1),size(img,2)).';
            rgbaVector=uint32(rgbaImg(:));
            
            fullScreenBufferedImage.setRGB(position(2),position(1),size(img,2),size(img,1),rgbaVector,0,size(img,2));
            sm.fullScreenWindows(id).icon.setImage(fullScreenBufferedImage);
            sm.fullScreenWindows(id).frame.pack();
            sm.fullScreenWindows(id).frame.repaint();
            sm.fullScreenWindows(id).frame.show();
        end
    end
    methods (Access=protected)
        function sm=openSingle(sm,id)
            if (~sm.isOpen(id))
                availableScreens=sm.screens;
                if (id<=length(availableScreens))
                    screenSize=availableScreens(id).size;
                    sm.fullScreenWindows(id).size=screenSize;

                    ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
                    screenDevices = ge.getScreenDevices();
                    screenDevice=screenDevices(id);

                    sm.fullScreenWindows(id).frame = javax.swing.JWindow(screenDevice.getDefaultConfiguration());
                    bounds = sm.fullScreenWindows(id).frame.getBounds();
                    sm.fullScreenWindows(id).bufferedImage=java.awt.image.BufferedImage(screenSize(2),screenSize(1),java.awt.image.BufferedImage.TYPE_INT_RGB);
                    sm.fullScreenWindows(id).icon = javax.swing.ImageIcon(sm.fullScreenWindows(id).bufferedImage);
                    label = javax.swing.JLabel(sm.fullScreenWindows(id).icon);
                    sm.fullScreenWindows(id).frame.getContentPane.add(label);
                    screenDevice.setFullScreenWindow(sm.fullScreenWindows(id).frame);
                    sm.fullScreenWindows(id).frame.setLocation(bounds.x, bounds.y);
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
                    sm.fullScreenWindows(id).frame.dispose();
                    sm.fullScreenWindows(id).frame=[];
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
            TF=id<=length(sm.fullScreenWindows) && ~isempty(sm.fullScreenWindows(id)) && ~isempty(sm.fullScreenWindows(id).frame);
        end
    end
end