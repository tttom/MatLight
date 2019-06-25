classdef QuadrantPhotoDiode < handle
    %  Quadrant Photo Diode class
    %  
    properties
    	meanIntensity=1;
	intensityGradient=[1 1 1];
    end
    methods
        function qpd=QuadrantPhotoDiode()
        	% some initialization here
	end
        function values=acquire(qpd)
            value=-1*ones(2);
        end
        function x=acquireX(qpd)
            x=diff(sum(aquire(qpd))/intensityGradient(1);
        end
        function y=acquireY(qpd)
            y=diff(sum(aquire(qpd).')/intensityGradient(2);
        end
        function z=acquireZ(qpd)
            z=(sum(sum(aquire(qpd))-qpd.meanIntensity)/qpd.intensityGradient(3);
        end
    end
    
end