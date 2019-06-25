classdef LightSource < handle
    % LightSource super class
    %
    properties (Abstract)
        targetPower;
    end
    
    methods
        function lightSource=LightSource()
        end
        function setWavelengths(source,wavelengths,relativePowers,gains)
        end
    end
end