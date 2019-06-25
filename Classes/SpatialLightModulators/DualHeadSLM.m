classdef DualHeadSLM < SLM
    % Dual Head SLM class
    %
    % constructor DualHeadSLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
    %
    % Example:
    % slm=DualHeadSLM(2,[1 1]/20,[128 128 512 512])
    %
    
    properties
        phaseChangeDueToAmplitudeModulation;
    end
    methods
        function slm=DualHeadSLM(displayNumber,referenceDeflectionFrequency,regionOfInterest)
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
            
            slm.twoPiEquivalent=1.0; % of dynamic range
            slm.stabilizationTime=0.30; % Conservative value
            slm.phaseChangeDueToAmplitudeModulation=0; %0.63/2;
        end
        function set.phaseChangeDueToAmplitudeModulation(slm,newPhaseChangeDueToAmplitudeModulation)
            slm.phaseChangeDueToAmplitudeModulation=newPhaseChangeDueToAmplitudeModulation;
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
            
            % Do some checking and clip values if needed
            if (max(abs(amplitudeValues(:)))>1)
                if (max(abs(amplitudeValues(:)))-1 > 4*eps(1))
                    %If not a rounding error, signal this.
                    logMessage('Warning: amplitude above 1, clipping it!');
                end
                amplitudeValues=min(amplitudeValues,1); %Convert to pixel value [0 1]
            end
            
            % Do any user defined action before we modulate the SLM
            slm.executeModulatePreFunctor(phaseValues,amplitudeValues);
            
            %Adjust for the correction function and reference
            phaseValues=mod((phaseValues + slm.referenceDeflectionAngle+slm.correctionFunctionAngle)./(2*pi)+0.5...
                -slm.phaseChangeDueToAmplitudeModulation*amplitudeValues,1); %Convert to pixel value [0 1)
            amplitudeValues=amplitudeValues.*slm.referenceDeflectionAbs.*slm.correctionFunctionAbs;
            
            %make sure that the amplitude is in the dynamic range of the SLM
            if (max(abs(amplitudeValues(:)))>1)
                amplitudeValues=amplitudeValues/max(eps(1),max(amplitudeValues(:))); %Convert to pixel value [0 1)
            end
            
            if (isscalar(amplitudeValues))
                amplitudeValues=amplitudeValues*ones(size(phaseValues));
            end
            imageForSLM=cat(3,zeros(size(phaseValues)),phaseValues,amplitudeValues);
            if (slm.twoPiEquivalent<=1)
                imageForSLM=slm.twoPiEquivalent*imageForSLM;
            else
                imageForSLM=min(slm.twoPiEquivalent*imageForSLM,1);
            end

            showImage(imageForSLM,[],slm.regionOfInterest(2),slm.regionOfInterest(1),slm.displayNumber);

            if (slm.stabilizationTime>0)
                pause(slm.stabilizationTime);%Make sure to wait long enough for the image to display on the SLM and stabilize
            end
            
            % Do any user defined action after we modulated the SLM
            slm.executeModulatePostFunctor(phaseValues,amplitudeValues,imageForSLM);
        end
        function delete(slm)
            if (~isempty(slm.displayNumber) && ishandle(slm.displayNumber) && strcmp(get(slm.displayNumber,'Type'),'axes'))
                % Close figure window
                close(get(slm.displayNumber,'Parent'));
            else
                % Close full screen window
                if (DefaultScreenManager.instance.screens(slm.displayNumber).open)
                    DefaultScreenManager.instance.close(slm.displayNumber);
                end
            end
        end
    end    
end

