%
% Singleton class that implements the ScreenManager abstract class.
% DefaultScreenManager provides an access point to the default full screen
% manager. If Cogent Graphics is installed, it will try to use that,
% otherwise it will rely on the built in Java full screen functionality.
%
% Properties:
%      screens: a vector with a struct for every detected screen at its index id.
%           struct fields:
%                id: the integer id (coincides with the index in the vector)
%                size: two-element vector indicating the height and width, respectively.
%                frameRate: the refresh rate of the screen in Hz.
%                open: boolean, true if the screen is currently open.
%                main: boolean, true if this screen holds the main desktop
%                description: English description of the screen.
%      openScreens: a list of id's of screens that are currently open.
%      mainScreen: the id of the main desktop screen.
%
% Methods:
%      instance(): returns the singleton instance of this class.
%      open(ids): opens the screens with id in the vector ids.
%      close(ids): closes the screens with id in the vector ids.
%      display(id,img,position): Opens the screen with id equal to the
%           first argument if it is not open yet. The image contained in
%           matrix img is displayed at the coordinates specified in
%           position. The image should be a 2D or 3D matrix. If it is a 3D
%           matrix, the third dimension should be 3, corresponding to the
%           red, green, and blue channel, respectively. If it is a 2D
%           matrix, the color is assumed to be gray. The position argument
%           contains the row and column position, followed by the height
%           and width of the region of interest on the screen. If position
%           is unspecified or empty, the top left corner of the screen is
%           assumed.
%
% Example usage:
%     sm=DefaultScreenManager.instance();
%     sm.screens.description % list all screen descriptions
%     sm.open(2); % Open screen 2
%     sm.display(2,getTestImage('boats'));
%     sm.display(2,getTestImage('boats'),[100 200 512 512]);
%     sm.display(sm.mainScreen,getTestImage('boats')); pause(2); sm.close(sm.mainScreen);
%     sm.close('all'); % or sm.delete();
%
% Note, you can always use the following to close all full-screens:
%     DefaultScreenManager.instance().delete();
%
%
classdef DefaultScreenManager < ScreenManager
    properties (SetAccess=private, Dependent=true)
        screens;
    end
    properties (Access=private)
        screenManager;
    end
    methods(Access=private)
        function sm=DefaultScreenManager()
            cogentGraphics=exist('cgloadlib','file');
            if (cogentGraphics)
                sm.screenManager=CogentGraphicsScreenManager.instance;
            else
                sm.screenManager=JavaScreenManager.instance;
            end
        end
    end
    methods (Static)
        function obj = instance()
            persistent uniqueInstance;
            if isempty(uniqueInstance)
                obj = DefaultScreenManager();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    methods
        function scrns=get.screens(sm)
            scrns=sm.screenManager.screens();
        end
        function sm=display(sm,id,img,position)
            switch(nargin)
                case 1
                    sm=sm.screenManager.display();
                case 2
                    sm=sm.screenManager.display(id);
                case 3
                    sm=sm.screenManager.display(id,img);
                otherwise
                    sm=sm.screenManager.display(id,img,position);
            end
        end
    end
    methods (Access=protected)
        function sm=openSingle(sm,id)
            sm=sm.screenManager.openSingle(id);
        end
        function sm=closeSingle(sm,id)
            sm=sm.screenManager.closeSingle(id);
        end
    end
end