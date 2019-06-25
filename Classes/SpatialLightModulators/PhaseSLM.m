classdef PhaseSLM < SLM
    % Phase SLM class
    %
    % constructor PhaseSLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
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
    % slm = PhaseSLM(2,[1 1]/20,[128 128 512 512])
    %
    methods
        function slm=PhaseSLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
            if (nargin<1)
                displayNumber=[];
            end
            if (nargin<2)
                referenceDeflectionFrequency=[];
            end
            if (nargin<3)
                regionOfInterest=[];
            end
            slm=slm@SLM(displayNumber,referenceDeflectionFrequency,regionOfInterest);
            
            slm.twoPiEquivalent=1.0; %Of dynamic range
            slm.stabilizationTime=0.30; % Conservative value
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
                phaseValues=complexFunction;
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
                if (abs(deflectionAngle)>pi/8),
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
            
            imageForSLM=cat(3,imageForSLM,imageForSLM,ones(size(imageForSLM)));%Set blue channel to 1 so that it can be used on a dual head setup as well
%             imageForSLM=repmat(imageForSLM,[1 1 3]);
            showImage(imageForSLM,[],slm.regionOfInterest(2),slm.regionOfInterest(1),slm.displayNumber);

            if (slm.stabilizationTime>0)
                pause(slm.stabilizationTime);%Make sure to wait long enough for the image to display on the SLM and stabilize
            end
            
            % Do any user defined action after we modulated the SLM
            slm.executeModulatePostFunctor(phaseValues,amplitudeValues,imageForSLM);
        end
        function delete(slm)
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
