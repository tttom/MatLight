%
% Abstract singleton class for full screen managers.
%
% Implementing classes: JavaScreenManager, CogentGraphicsScreenManager, and
% DefaultScreenManager.
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
%     sms.screens.description % list all screen descriptions
%     sm.open(2); % Open screen 2
%     sm.display(2,getTestImage('boats'));
%     sm.display(2,getTestImage('boats'),[100 200 512 512]);
%     sm.display(sm.mainScreen,getTestImage('boats')); pause(2); sm.close(sm.mainScreen);
%     sm.close('all'); % or sm.delete();
%
classdef ScreenManager < Singleton
    properties (SetAccess=private, Dependent=true)
        mainScreen;
        openScreens;
    end
    properties (Abstract, SetAccess=private, Dependent=true)
        screens;
    end
    methods (Access=protected)
        function sm=ScreenManager()
        end
    end
    methods
        function sm=delete(sm)
            sm.close('all');
        end
        function sm=open(sm,ids)
            if (isempty(ids) || (ischar(ids) && strcmpi(ids,'all')))
                ids=sm.openScreens.id;
            end
            for (id=ids)
                sm=sm.openSingle(id);
            end
        end
        function sm=close(sm,ids)
            if (nargin<2 || isempty(ids) || (ischar(ids) && strcmpi(ids,'all')))
                ids=sm.openScreens;
            end
            for (id=ids)
                sm=sm.closeSingle(id);
            end
        end
        function openScrns=get.openScreens(sm)
            openScrns=cell2mat({sm.screens(cell2mat({sm.screens.open})==true).id});
        end
        function scrn=get.mainScreen(sm)
            scrn = find(cell2mat({sm.screens.main})==true);
        end
    end
    methods (Abstract)
        sm=display(sm,id,img,position);
    end
    methods (Abstract, Access=protected)
        sm=openSingle(sm,id);
        sm=closeSingle(sm,id);
    end
end