classdef BNSPhaseSLM < SLM
    % BNS Phase SLM class
    % Specially for the Bolder Linear System (BNS) PCI Phase SLM 
    %
    % constructor BNSPhaseSLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
    %
    % displayNumber: either an integer indicating the output port of the
    %      graphics card, or a reference to a Matlab axes
    % referenceDeflectionFrequency: vector of two elements indicating the
    %      vertical and horizontal grating frequency in units of pixels^-1,
    %      a third vector element can indicate additional defocus
    % regionOfInterest: is the part of the SLM to be used specified as [down right height width]
    %
    % All arguments are optional.
    %
    % Example
    % slm=BNSPhaseSLM(-1,[1 1 0]/20,[0 0 512 512])
    %                   The first param must be -1 for PCI BNS SLM;
    %                   Ortherwise it is for DVI BNS SLM;
    % By Mingzhou Chen @ 07/02/2014 mc225@st-andrews.ac.uk Ver. 1.0
    
    properties
        currentBoard; %the board need to be controlled;
    end
    properties (SetAccess = protected)
        boardNumber; %The number of boards connected to the PC
        referenceDeflectionExt = [];
        isWavefrontCorrectionLoaded;
        slmX;
        slmY;
    end
    
    methods
        function slm=BNSPhaseSLM(displayNumber,dlldir,LUTFileName,referenceDeflectionFrequency)
            if nargin<1
                displayNumber = -1; %PCI BNS SLM;
            end
            
            if (nargin<3)
                dlldir = 'C:\PCIe16MatlabSDK\';           %location of all BNS SDK files;
                LUTFileName = 'C:\PCIe16MatlabSDK\LUT_Files\slm3876_at1064_P16.lut'; %need to copy the lut file into this directory;
            end
            
            if (nargin<4)
                referenceDeflectionFrequency = [];
            end
            
            slm=slm@SLM(displayNumber,referenceDeflectionFrequency,[0 0 512 512]);  %-1 means slm is not as a monitor;
            if libisloaded('Interface')
                unloadlibrary('Interface');
            end
            loadlibrary([dlldir 'Interface.dll'], [dlldir 'Interface.h']);
            slm.boardNumber = calllib('Interface', 'Constructor', 1); %constructor
            if slm.boardNumber ==0
                disp('No BNS board can be found, please check connections and retry!');
            else
                if slm.boardNumber>=1
                    slm.currentBoard = 1;
                    calllib('Interface', 'SetTrueFrames', 3);%setup tru
                    calllib('Interface', 'SLMPower', 1); %power on SLM;
                    calllib('Interface', 'LoadLUTFile', slm.currentBoard, LUTFileName);
                end
            end         
            
            %reset slm roi and maxSize; have to change here for the SLM with other resolution;
            slm.slmX = slm.maxSize(1);
            slm.slmY = slm.maxSize(2);
            
            slm.isWavefrontCorrectionLoaded = 0; 
            slm.twoPiEquivalent=1.0; %Of dynamic range
            slm.stabilizationTime=0.1; % Conservative value      
            
            slm=updateReferenceDeflectionExt(slm);
        end
        
       function boardNumber=get.boardNumber(slm)
            boardNumber=slm.boardNumber;
       end
       function slm=set.boardNumber(slm,boardNumber)
            slm.boardNumber = boardNumber;
       end
       
       function currentBoard=get.currentBoard(slm)
            currentBoard=slm.currentBoard;
       end
       function slm=set.currentBoard(slm,currentBoard)
           if currentBoard <=0
               slm.currentBoard = 1;
           elseif currentBoard <= slm.boardNumber 
               slm.currentBoard = currentBoard;
           else               
               slm.currentBoard = 1;
           end
       end
       
       function modulate(slm,complexPupilFunction,amplitudePupilFunction)
            if (nargin<2 || isempty(complexPupilFunction))
                complexPupilFunction=0;
            else
                if (isa(complexPupilFunction,'function_handle'))
                    complexPupilFunction=complexPupilFunction(slm.X,slm.Y);
                end
            end
            if (nargin<3 || isempty(amplitudePupilFunction))
                phaseValues=angle(complexPupilFunction);
                amplitudeValues=abs(complexPupilFunction); %Assume that the reference deflection is just that, no amplitude mod allowed there.
            else
                phaseValues=complexPupilFunction;
                %Handle phase and amplitude separately
                if (isa(amplitudePupilFunction,'function_handle'))
                    amplitudeValues=amplitudePupilFunction(slm.X,slm.Y);
                else
                    amplitudeValues=amplitudePupilFunction;
                end
            end
            
            % Do any user defined action before we modulate the SLM
            slm.executeModulatePreFunctor(phaseValues,amplitudeValues);
            
            %Adjust for the correction function and reference
            phaseValues=mod((phaseValues + slm.referenceDeflectionAngle+slm.correctionFunctionAngle)./(2*pi)+0.5,1); %Convert to pixel value [0 1)
          %  phaseValues=((phaseValues + slm.referenceDeflectionAngle+slm.correctionFunctionAngle)./(2*pi)+0.5); %Convert to pixel value [0 1)
            amplitudeValues=amplitudeValues.*slm.referenceDeflectionAbs.*slm.correctionFunctionAbs;
            
            %make sure that the amplitude is in the dynamic range of the SLM
            if (max(abs(amplitudeValues(:)))>1)
                if (max(abs(amplitudeValues(:)))-1 > 4*eps(1))
                    %If not a rounding error, signal this.
                    logMessage('Warning: amplitude above 1, clipping it!');
                end
                amplitudeValues=min(amplitudeValues,1); %Convert to pixel value [0 1]
            end
            
            %Modulate the amplitude with the phase SLM
            phaseDeviationFromMeanForAmplitude=acos(amplitudeValues)./(2*pi); %this is still the absolute value
                        
            %Select pixels for binary coding
            if (~all(slm.referenceDeflectionFrequency(1:2)==0))
                deflectionAngle=mod(atan2(slm.referenceDeflectionFrequency(2),slm.referenceDeflectionFrequency(1))+pi/4,pi/2)-pi/4;
                if (abs(deflectionAngle)<=pi/8)
                    % Dump intensity in the highest diagonal orders
                    selectedPixels=1.0*(mod(slm.X+slm.Y,2)>0);
                else
                    % Don't dump intensity diagonally because the first order
                    % deflection is too close to the diagonal
                    if (abs(slm.referenceDeflectionFrequency(1))>abs(slm.referenceDeflectionFrequency(2)))
                        selectedPixels=1.0*(mod(slm.X,2)>0);
                    else
                        selectedPixels=1.0*(mod(slm.Y,2)>0);
                    end
                end
            else
                % We are doing zeroth order modulation
                selectedPixels=1.0*(mod(slm.X+slm.Y,2)>0); %dump intensity at highest orders
%                 selectedPixels=mod(slm.X/slm.dumpRelativeSpatialFrequency(1)+slm.Y/slm.dumpRelativeSpatialFrequency(2),1);
            end
            
            %Convert the absolute value to a signed value. Half of the pixels should be biased upwards and the other half downwards
            phaseDeviationFromMeanForAmplitude=phaseDeviationFromMeanForAmplitude.*(2*selectedPixels-1);
            
            imageForSLM=mod(phaseValues+phaseDeviationFromMeanForAmplitude,1);
            if (slm.twoPiEquivalent<=1)
                imageForSLM=slm.twoPiEquivalent*imageForSLM;
            else
                imageForSLM=max(min(slm.twoPiEquivalent*(imageForSLM-.5)+.5,1),0); %Center around midpoint of the dynamic range of the SLM
            end
            
%             imageForSLM=cat(3,imageForSLM,imageForSLM,ones(size(imageForSLM)));%Set blue channel to 1 so that it can be used on a dual head setup as well
% %             imageForSLM=repmat(imageForSLM,[1 1 3]);
%             showImage(imageForSLM,[],slm.regionOfInterest(2),slm.regionOfInterest(1),slm.displayNumber);

            normmask = uint16(imageForSLM*(2^16)); %turn into a 16bit array; 
             %typecast will arrange one pixel as [low8bit high8bit];
            slmmask= libpointer('uint8Ptr',typecast(reshape(normmask.',1,[]),'uint8'));  %need this .'
            
            calllib('Interface', 'WriteImage', slm.currentBoard, slmmask);

            
            if (slm.stabilizationTime>0)
                pause(slm.stabilizationTime);%Make sure to wait long enough for the image to display on the SLM and stabilize
            end
            
            % Do any user defined action after we modulated the SLM
            slm.executeModulatePostFunctor(phaseValues,amplitudeValues,imageForSLM);
        end
        
        function modulate1(slm,mask,isWavefrontCorr)
            isMaskOk = 1;
            if nargin < 2
                mask = zeros(slm.slmX,slm.slmY);
            end
            if nargin < 3
                isWavefrontCorr = 0;
            end
            if ischar(mask)
                if exist(mask,'file')==2 %load from tiff file;
                    slmmask = libpointer('uint8Ptr', zeros(slm.slmX*slm.slmY*2,1));
                    calllib('Interface', 'ReadTIFF', mask, slmmask);
                else
                    fprintf('Can not find the mask tiff file:\n %s.\n',mask);
                    isMaskOk = 0;
                end
            else %load from current work space, which is prepared, 512x512 image;
                if mask == 0 ; %block SLM
                    mask = zeros(slm.slmX,slm.slmY);
                elseif mask == 1; %All SLM on;
                    mask = ones(slm.slmX,slm.slmY);
                end
                
                if length(mask(:))~=slm.slmX*slm.slmY
                    logMessage(sprintf('Warning: Mask must be %dx%d pixels!',slm.slmX, slm.slmY));
                    isMaskOk = 0;
                end
                
                if isMaskOk
                    if ~isreal(mask)
                        disp('Amplitude will be encoded in phase on SLM!');
                        mask = angle(mask.*slm.referenceDeflection+(1-abs(mask)).*slm.referenceDeflectionExt);
                    else
                        mask = angle(exp(1i*mask).*slm.referenceDeflection);
                    end
                    
                    normmask = (mask+pi)./2/pi; %normalized by 2pi;
                    normmask = uint16(normmask*(2^16)); %turn into a 16bit array;
                    %typecast will arrange one pixel as [low8bit high8bit];
                    slmmask= libpointer('uint8Ptr',typecast(reshape(normmask.',1,[]),'uint8'));  %need this .'
                end
            end
            
            if isMaskOk            
                if isWavefrontCorr %wavefront correction;
                    calllib('Interface', 'WriteCal', slm.currentBoard, slmmask);
                    logMessage('Warning: New wavefront correction has been applied...');
                    slm.isWavefrontCorrectionLoaded = 1;
                else
                    calllib('Interface', 'WriteImage', slm.currentBoard, slmmask);
                end
            end               
            
            if (slm.stabilizationTime>0)
                pause(slm.stabilizationTime);%Make sure to wait long enough for the image to display on the SLM and stabilize
            end
        end
        function slm=updateReferenceDeflectionExt(slm)
            refFreq = slm.referenceDeflectionFrequency;
            slm.referenceDeflectionFrequency = [refFreq(1) -refFreq(2) refFreq(3)];
            slm.referenceDeflectionExt=slm.referenceDeflection;
            slm.referenceDeflectionFrequency = refFreq;
        end
        function delete(slm)
            if ~libisloaded('Interface')
                logMessage('No SLM dll has been loaded!');
            else
                calllib('Interface', 'SLMPower', 0); %power off SLM;
                calllib('Interface', 'Deconstructor'); %desctructor;
                unloadlibrary('Interface');
                logMessage('SLM has been deleted!');
            end
            if (~isempty(slm.displayNumber) && ishandle(slm.displayNumber) && strcmp(get(slm.displayNumber,'Type'),'axes'))
                close(get(slm.displayNumber,'Parent'));
            else
                if (DefaultScreenManager.instance.screens(slm.displayNumber).open)
                    DefaultScreenManager.instance.close(slm.displayNumber);
                end
            end
        end
    end    
end
