classdef SLM < handle
    % Abstract SLM class
    %
    % constructor SLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
    %
    % The default displayNumber is selected by [], which is the first screen that is not the currently used screen, or if there is only one, a pop-up window.
    %
    % Example
    % slm=PhaseSLM(2,[1 1]/20,[128 128 512 512])
    % slm=DualHeadSLM(2,[1 1]/20,[128 128 512 512])
    %
    properties
        referenceDeflectionFrequency=[1 1 0]./20; % [y, x, w20], the frequency in cycles/pixel, and the defocus in wavelengths per pixel^2
        twoPiEquivalent=1.0; %Of dynamic range
        stabilizationTime=0.0;
        
        modulatePreFunctor=[];
        modulatePostFunctor=[];
    end
    properties(Dependent = true)
        regionOfInterest; % [toprow leftcolumn rows cols]
        correctionFunction=[]; %Complex matrix of the size of the region of interest
    end
    properties (SetAccess = protected)
        maxSize=[1 1]; % [rows cols]
        displayNumber;
        referenceDeflection=[]; %Complex matrix of the size of the region of interest
        referenceDeflectionAngle=[]; %Real matrix with the deflection phase in [-pi, pi)
        referenceDeflectionAbs=[]; %The amplitude of the reference deflection. Can be a matrix or a scalar to represent a constant matrix.
    end
    properties (Access = protected)
        regionOfInterestOffset=[0 0];
        X;
        Y;
        
        correctionFunctionComplex;
        correctionFunctionAngle;
        correctionFunctionAbs;
    end
    
    methods
        function slm=SLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
            if (nargin<1 || isempty(displayNumber))
                % If the display isn't specified, pick something else than the main display
                ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
                screenDevices = ge.getScreenDevices();
                mainScreen = ge.getDefaultScreenDevice();
                displayNumber=1;
                while (screenDevices(displayNumber)==mainScreen && displayNumber<length(screenDevices))
                    displayNumber=displayNumber+1;
                end
                if (screenDevices(displayNumber)==mainScreen)
                    figH=figure();
                    displayNumber=axes('Parent',figH);
                    image(zeros(600,800),'Parent',displayNumber);
                end
            end
            slm.displayNumber=displayNumber;
            slm.maxSize=slm.detectResolution(displayNumber);
            slm.correctionFunction=1;
            if (~isempty(slm.maxSize)),
                if (nargin<2 || isempty(referenceDeflectionFrequency))
                    referenceDeflectionFrequency=[0 0 0];
                end
                if (nargin<3 || isempty(regionOfInterest))
                    regionOfInterest=[0 0 slm.maxSize];
                end
                slm.regionOfInterest=regionOfInterest;
                if (length(referenceDeflectionFrequency)<3)
                    referenceDeflectionFrequency(3)=0;
                end
                slm.referenceDeflectionFrequency=referenceDeflectionFrequency(1:3);
            else
               error('Display %u not found!',displayNumber);
            end
        end
        function slm=set.regionOfInterest(slm,newROI)
            if (length(newROI)<3)
                newROI(3:4)=slm.maxSize;
            end
            %Clip the region of interest to the active area of the SLM
            newROI(1:2)=min(max(0,newROI(1:2)),slm.maxSize-1);
            newROI(3:4)=min(newROI(1:2)+max(2,newROI(3:4)),slm.maxSize)-newROI(1:2);
            
            slm.regionOfInterestOffset=newROI(1:2);
            slm=updateReferenceDeflection(slm,newROI(3:4));
        end
        function roi=get.regionOfInterest(slm)
            roi=[slm.regionOfInterestOffset size(slm.referenceDeflection)];
        end
        function slm=set.referenceDeflectionFrequency(slm,newReferenceDeflectionFrequency)
            if (any(abs(newReferenceDeflectionFrequency(1:2))>=1/2)),
                logMessage('The reference deflection is (%d,%d) [cycles/slm-pixel], while it should generally be smaller than 1/2th in absolute value',newReferenceDeflectionFrequency(1:2));
            end
            if (length(newReferenceDeflectionFrequency)<3)
                newReferenceDeflectionFrequency(3)=0;
            end
            slm.referenceDeflectionFrequency=newReferenceDeflectionFrequency;
            slm=updateReferenceDeflection(slm);
        end
        function slm=set.correctionFunction(slm,newCorrectionFunction)
            if (isa(newCorrectionFunction,'function_handle'))
            	newCorrectionFunction=newCorrectionFunction(slm.X,slm.Y);
            end
            if (any(size(newCorrectionFunction)>1) && any(size(newCorrectionFunction)~=slm.regionOfInterest(3:4)))
                newROI=[floor((1+slm.maxSize-size(newCorrectionFunction))/2) size(newCorrectionFunction)];
                if (all(newROI(1:2)>=0))
                    logMessage('Updating region of interest of SLM because the specified aberration correction has a different size.');
                    slm.regionOfInterest=newROI;
                else
                    error('The specified aberration is larger than the SLM.');
                end
            end
            
            slm.correctionFunctionComplex=newCorrectionFunction;
            slm.correctionFunctionAngle=angle(slm.correctionFunction);
            slm.correctionFunctionAbs=abs(slm.correctionFunction);
        end
        function correctionFunction=get.correctionFunction(slm)
            correctionFunction=slm.correctionFunctionComplex;
        end
    end
    methods (Abstract)
        modulate(slm,complexPupilFunction,amplitudePupilFunction);
    end
    methods (Access = protected)
        function slm=executeModulatePreFunctor(slm,phaseValues,amplitudeValues)
            if (~isempty(slm.modulatePreFunctor))
                inputArgs = {};
                if (nargin(slm.modulatePreFunctor)>=1),
                    inputArgs{1} = amplitudeValues.*exp(1i*phaseValues); % complexPupilFunction
                end
                if (nargout(slm.modulatePreFunctor)<=0)
                    slm.modulatePreFunctor(inputArgs{:});
                else
                    slm = slm.modulatePreFunctor(inputArgs{:});
                end
                end
            end
        function slm=executeModulatePostFunctor(slm,phaseValues,amplitudeValues,imageOnSLM)
            if (~isempty(slm.modulatePostFunctor))
                inputArgs = {};
                if (nargin(slm.modulatePostFunctor)>=1),
                    inputArgs{1} = amplitudeValues.*exp(1i*phaseValues); % complexPupilFunction
                end
                if (nargin(slm.modulatePostFunctor)>=2),
                    inputArgs{2} = imageOnSLM;
                end
                if (nargout(slm.modulatePostFunctor)<=0),
                    slm.modulatePostFunctor(inputArgs{:});
                else
                    slm = slm.modulatePostFunctor(inputArgs{:});
                end
            end
        end
    end
    methods (Access = private)
        function slm=updateReferenceDeflection(slm,regionOfInterestSize)
            if (nargin<2 || isempty(regionOfInterestSize))
                regionOfInterestSize=size(slm.referenceDeflection);
            end
            regionOfInterestCenter=floor(regionOfInterestSize./2)+1;
            [slm.X,slm.Y]=meshgrid([1:regionOfInterestSize(2)]-regionOfInterestCenter(2),[1:regionOfInterestSize(1)]-regionOfInterestCenter(1));
            slm.referenceDeflectionAngle=2*pi*(slm.X*slm.referenceDeflectionFrequency(2)+slm.Y*slm.referenceDeflectionFrequency(1)+slm.referenceDeflectionFrequency(3)*(slm.X.^2+slm.Y.^2));
            slm.referenceDeflection=exp(1i*slm.referenceDeflectionAngle);
            slm.referenceDeflectionAbs=1.0;
        end
    end
    methods (Static = true)
        function displaySize=detectResolution(displayNumber)
            ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
            screenDevices = ge.getScreenDevices();
            if (isnumeric(displayNumber) && (displayNumber-round(displayNumber)<eps('single')) && displayNumber>0 && displayNumber<=length(screenDevices)),
                displaySize(1) = screenDevices(displayNumber).getDisplayMode().getHeight();
                displaySize(2) = screenDevices(displayNumber).getDisplayMode().getWidth();
            else
                if (ishandle(displayNumber) && strcmpi(get(displayNumber,'Type'),'axes'))
                    img=findobj(displayNumber,'Type','image');
                    if isempty(img)
                        img=image(zeros(300,400),'Parent',displayNumber);
                    end
                    displaySize=size(get(img,'CData'));
                    displaySize=displaySize(1:2);
                else
                    displaySize=[];
                end
            end
        end
    end
    
end

