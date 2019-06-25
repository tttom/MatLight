% A GUI to use speckle as a wavelength meter.
% Load a mat file created by analyzeSpeckleImages, and select the camera in use from the menu.
%
function waveMeter(calibrationFileName)  
    % Find the screen size and put the application window in the center
    set(0,'units','pixels');
    screenSize=get(0,'screensize');
    screenSize=screenSize([3 4]);
    windowSize=[600 145];
    fig=figure('Units','pixels','Position',[floor(screenSize/2)-floor(windowSize/2) windowSize],'MenuBar','none','Color',[0 0 0],'NumberTitle','off','CloseRequestFcn',@exit);
    textDisplay=uicontrol('Parent',fig,'Style','text','Position',[25 25 550 95],'FontName','Arial','FontSize',64,'FontWeight','bold','BackgroundColor',[0 0 0]);
    setUserData('textDisplay',textDisplay);
    menuFile=uimenu(fig,'Label','File');
    menuLoadCalibration=uimenu(menuFile,'Label','Load Calibration...','Callback',@selectAndLoadCalibration);
    menuExit=uimenu(menuFile,'Label','Exit','Callback',@exit);
    menuCam=uimenu(fig,'Label','Camera');
    detectedCameras=detectCameras();
    setUserData('detectedCameras',detectedCameras);
    for (camIdx=1:length(detectedCameras))
        uimenu(menuCam,'Label',detectedCameras(camIdx).description,'Callback',@(obj,event) updateCam(detectedCameras(camIdx)));
    end
    menuSignal=uimenu(fig,'Label','Signal','Callback',@openSignalWindow);
    % Window intialized now
    
    loadStatus();
    updateStatus('version',0.10);
    % Initial values
    setUserData('wavelengths',[]);
    setUserData('regionOfInterest',[]);
    setUserData('isClosing',false);
    setUserData('imageAxes',[]);
    setUserData('cam',[]);
    
    if (nargin<1 || isempty(calibrationFileName))
        calibrationFileName=getStatus('calibrationFileName',[]);
    end
    updateStatus('calibrationFileName',calibrationFileName);
    
    % Hook up a camera
    selectedCamera={};
    selectedCamera.type=getStatus('selectedCameraType',[]);
    selectedCamera.index=getStatus('selectedCameraIndex',[]);
    updateCam(selectedCamera);
    
    % Load the calibration file or ask the user
    loadCalibration(calibrationFileName);
    
    % Loop until stopped
    while (~getUserData('isClosing')),
        if (~isempty(getUserData('cam')))
            updateDisplay();
            drawnow();
        else
            pause(.020);
        end
    end
    
    % Close any extra windows if needed
    closeSignalWindow();
    
    %Shut the main window
    closereq();
end

function setFigureTitle()
    wavelengths=getUserData('wavelengths');
    if (~isempty(wavelengths))
        precision=max(diff(wavelengths));
        significantDigits=min(4,max(0,-log10(precision)-9));
        significantDigitsPrecision=1;
        titleString=sprintf(sprintf('%%0.%0.0ff nm to %%0.%0.0ff nm, %%0.%0.0ff pm resolution',[[1 1]*significantDigits significantDigitsPrecision]),1e9*[min(wavelengths) max(wavelengths) 1e3*precision]);
    else
        titleString='no calibration data loaded';
    end
    set(getBaseFigureHandle(),'NumberTitle','off','Name',strcat('waveMeter - ',titleString));
end

function updateDisplay()
    textDisplay=getUserData('textDisplay');
    
    [currentWavelength status]=determineWavelength();
    
    switch (status)
        case 'signalTooLow'
            set(textDisplay,'ForegroundColor',[1 0 0],'String','LOW SIGNAL');
        case 'signalTooHigh'
            set(textDisplay,'ForegroundColor',[.2 0.2 1],'String','TOO BRIGHT');
        case 'outOfRange'
            set(textDisplay,'ForegroundColor',[.4 0.8 0.4],'String','OUT OF RANGE');
        case 'noCalibrationLoaded'
            set(textDisplay,'ForegroundColor',[0 0.8 0],'String','NOT CALIB.');
        otherwise
            set(textDisplay,'ForegroundColor',[1 1 1],'String',sprintf('%0.4f nm',currentWavelength*1e9));
    end
end

function [wavelength status]=determineWavelength()
    shortestIntegrationTime=10e-6;
    longestIntegrationTime=.5;
    
    status=[]; %default
    wavelength=[]; %default
    
    % capture image
    cam=getUserData('cam');
    img=cam.acquire();
    
    while ((max(img(:))>=1 && cam.integrationTime>shortestIntegrationTime*2) || (max(img(:))<.5 && cam.integrationTime<longestIntegrationTime/2))
        if (max(img(:))>=1)
            cam.integrationTime=cam.integrationTime/2;
        else
            cam.integrationTime=cam.integrationTime*2;
        end
%         cam.defaultNumberOfFramesToAverage=max(1,floor(.03/cam.integrationTime));
%         logMessage('Integration time %f ms',cam.integrationTime*1e3);
        img=cam.acquire();
    end
    setUserData('cam',cam);
    
    %Display output if required
    imageAxes=getUserData('imageAxes');
    if (~isempty(imageAxes))
        if (ishandle(imageAxes) && strcmpi(get(imageAxes,'BeingDeleted'),'off'))
            showImage(img,-1,[],[],imageAxes);
        else
            setUserData('imageAxes',[]);
        end
    end
    
    % Check the input
    if (max(img(:))>=1)
        status='signalTooHigh';
    end
    if (max(img(:))<=.5)
        status='signalTooLow';
    end
    wavelengths=getUserData('wavelengths');
    if (isempty(wavelengths))
        status='noCalibrationLoaded';
    end
    
    if (isempty(status)) % no errors
        calibration=getUserData('calibration');
        
        wavelength = determineWavelengthFromSpeckleImage(img,calibration);
        %TODO Check match error and output 'outOfRange' if needed
        
        status='ok'; % no errors
    end
end

function loadCalibration(fullFileName)
    try
        calibration=load(fullFileName);
        setUserData('calibration',calibration);

        cam=getUserData('cam');
        if (~isempty(cam))
            try
                cam.regionOfInterest=calibration.regionOfInterest;
            catch Exc
                logMessage('Region of interest invalid, please use the same camera!');
                cam=[];
            end
            setUserData('cam',cam);
        end
        wavelengths=sort(unique(calibration.trainingWavelengths));
        setUserData('wavelengths',wavelengths);

        logMessage('Loaded calibration data for %d wavelengths from %0.1f nm to %0.1f nm',[length(wavelengths) 1e9*[min(wavelengths) max(wavelengths)]]);
        updateStatus('calibrationFileName',fullFileName);
        setFigureTitle();
    catch Exc
        selectAndLoadCalibration();
    end
    
    setFigureTitle();
end

function cameras=detectCameras()
    cameras=[];
    cameras(1).type='DummyCam';
    cameras(1).index=1;
    cameras(1).description='Dummy Cam';
    cameras(2).type='DummyCam';
    cameras(2).index=2;
    cameras(2).description='Dummy Cam HD';
    
    hwInfo=imaqhwinfo();
    if (any(strcmpi(hwInfo.InstalledAdaptors,'gige')))
    	hwInfoGigE=imaqhwinfo('gige');
        for (camIdx=1:length(hwInfoGigE.DeviceIDs))
            cameras(end+1).type='BaslerGigE';
            cameras(end).index=hwInfoGigE.DeviceIDs{camIdx};
            cameras(end).description='Basler GigE';
            if (length(hwInfoGigE.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoGigE.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'andor')))
    	hwInfoAndor=imaqhwinfo('andor');
        for (camIdx=1:length(hwInfoAndor.DeviceIDs))
            cameras(end+1).type='Andor';
            cameras(end).index=hwInfoAndor.DeviceIDs{camIdx};
            cameras(end).description='Andor';
            if (length(hwInfoGigE.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoAndor.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'winvideo')))
    	hwInfoIC=imaqhwinfo('winvideo');
        for (camIdx=1:length(hwInfoIC.DeviceIDs))
            cameras(end+1).type='ImagingSource';
            cameras(end).index=hwInfoIC.DeviceIDs{camIdx};
            cameras(end).description='Imaging Source';
            if (length(hwInfoIC.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoIC.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'avtmatlabadaptor64_r2009b')))
    	hwInfoPC=imaqhwinfo('avtmatlabadaptor64_r2009b');
        for (camIdx=1:length(hwInfoPC.DeviceIDs))
            cameras(end+1).type='PikeCam';
            cameras(end).index=hwInfoPC.DeviceIDs{camIdx};
            cameras(end).description='PikeCam';
            if (length(hwInfoPC.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoPC.DeviceIDs{camIdx}));
            end
        end
    end
    if (exist('vcapg2'))
        numberOfCards=vcapg2();
        for cardIdx=1:numberOfCards,
            cameras(end+1).type='DirectShow';
            cameras(end).index=cardIdx;
            cameras(end).description='Win Direct Show';
            if (numberOfCards>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',cardIdx));
            end
        end
    end
end

%
% Menu callbacks
%
function selectAndLoadCalibration(obj,event)
    [fileName,pathName]=uigetfile('*.mat','Select calibration file...');
    figure(getBaseFigureHandle()); % Bring the display back to the foreground
    if (~isempty(fileName) && ischar(fileName))
        fullFileName=strcat(pathName,'/',fileName);
        loadCalibration(fullFileName);
    else
        logMessage('No file selected, keeping current calibration data.');
    end
end
function exit(obj,event)
    setUserData('isClosing',true);
end
function openSignalWindow(obj,event)
    imageAxes=getUserData('imageAxes');
    if (isempty(imageAxes))
        cam=getUserData('cam');
        mainFigHandle=getBaseFigureHandle();
        windowOffset=get(mainFigHandle,'Position');
        windowOffset=windowOffset(1:2)+windowOffset(3:4)-[0 cam.regionOfInterest(3)];
        signalFig=figure('Name','Live Signal','Units','pixels','Position',[windowOffset cam.regionOfInterest([4 3])],'NumberTitle','off','UserData',struct('mainFigure',mainFigHandle));
        imageAxes=axes('Parent',signalFig,'Units','normalized','Position',[0 0 1 1]);
        setUserData('imageAxes',imageAxes);
    else
        %Already open, user probably want to close it instead
        closeSignalWindow();
    end
end
function closeSignalWindow()
    imageAxes=getUserData('imageAxes');
    if (~isempty(imageAxes))
        if (ishandle(imageAxes) && strcmpi(get(imageAxes,'BeingDeleted'),'off'))
            close(get(imageAxes,'Parent'));
        end
    end
end
function updateCam(selectedCamera)
    try
        switch(selectedCamera.type),
            case 'BaslerGigE'
                %Basler GigE
                cam=BaslerGigECam(selectedCamera.index);
            case 'Andor'
                %Andor
                logMessage('Andor camera not implemented');
            case 'ImagingSource'
                %Imaging Source
                cam=ImagingSourceCam(selectedCamera.index);
            case 'DirectShow'
                %Direct Show
                cam=DirectShowCam(selectedCamera.index);
            case 'PikeCam'
                %Pike Show
                cam=PikeCam(selectedCamera.index);
            otherwise
                %DummyCam
                switch (selectedCamera.index)
                    case 2
                        imgSize=[480 640];
                    otherwise
                        imgSize=[240 320];
                end
                cam=DummyCam(imgSize);
        end
        try
            cam.regionOfInterest=getUserData('regionOfInterest');
            setUserData('cam',cam);
            updateStatus('selectedCameraType',selectedCamera.type);
            updateStatus('selectedCameraIndex',selectedCamera.index);
        catch Exc
            logMessage('Region of interest invalid, please use the same camera!');
            cam=[];
        end
    catch Exc
        logMessage('Could not set camera!');
    end
end

%
% Data access functions
%
%
% User data function are 'global' variables linked to this GUI, though not
% persistent between program shutdowns.
%
function value=getUserData(fieldName)
    fig=getBaseFigureHandle();
    userData=get(fig,'UserData');
    value=userData.(fieldName);
end
function setUserData(fieldName,value)
    if (nargin<2)
        value=fieldName;
        fieldName=inputname(1);
    end
    fig=getBaseFigureHandle();
    userData=get(fig,'UserData');
    userData.(fieldName)=value;
    set(fig,'UserData',userData);
end
function fig=getBaseFigureHandle()
    fig=gcbf();
    if (isempty(fig))
        fig=gcf();
    end
    userData=get(fig,'UserData');
    if (isfield(userData,'mainFigure'))
        fig=userData.mainFigure;
    end
end
function loadStatus()
    fullPath=mfilename('fullpath');
    try
        load(strcat(fullPath,'.mat'),'status');
    catch Exc
    end
    if (~exist('status','var'))
        status={};
        status.version=-1;
    end
    
    setUserData('status',status);
end
function value=getStatus(fieldName,defaultValue)
    status=getUserData('status');
    if (isfield(status,fieldName))
        value=status.(fieldName);
    else
        if (nargin<2)
            defaultValue=[];
        end
        value=defaultValue;
    end
end
function updateStatus(fieldName,value)
    status=getUserData('status');
    status.(fieldName)=value;
    setUserData('status',status);
    
    saveStatus();
end
function saveStatus()
    status=getUserData('status');
    
    fullPath=mfilename('fullpath');
    save(strcat(fullPath,'.mat'),'status');
end

