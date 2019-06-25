% 
% A class representing a dimension, including units
%
classdef Dim
    properties (SetAccess=private)
        value;
        unit;
    end
    properties(Constant)
        y=1e-24; yocto=Metric.y;
        z=1e-21; zepto=Metric.z;
        a=1e-18; atto=Metric.a;
        f=1e-15; femto=Metric.f;
        p=1e-12; pico=Metric.p;
        n=1e-9; nano=Metric.n;
        u=1e-6; micro=Metric.u;
        m=1e-3; milli=Metric.m;
        
        c=1e-2; centi=Metric.c;
        d=1e-1; deci=Metric.d;
        da=1e1; deca=Metric.da;
        h=1e2; hecto=Metric.h;
        
        k=1e3; kilo=Metric.k;
        M=1e6; mega=Metric.M;
        G=1e9; giga=Metric.G;
        T=1e12; tera=Metric.T; 
        P=1e15; peta=Metric.P;
        E=1e18; exa=Metric.E;
        Z=1e21; zetta=Metric.Z;
        Y=1e24; yotta=Metric.Y;
    end
    properties(Constant,Access=private)
        prefixes1000={'y','z','a','f','p','n','u','m','','k','M','G','T','P','E','Z','Y'};
        prefixes10={'c','d','','da','h'};
    end
    methods
        function dim=Dim(dimensionStr,unitStr)
            if nargin<2 || isempty(unitStr),
                dim.value=regexp(dimensionStr,'%f');
            else
                dim.value=dimensionStr;
                dim.unit=parseUnit(unitStr);
            end
        end
        function str=disp(dim)
            exponent=floor(log10(dim.value));
            scaling=10^exponent;
            str=sprintf('%0.3f',dim.value/scaling)+prefix(scaling)+dim.unit;
        end
    end
    methods (Static=true)
        function prefixStr=prefix(value)
            exponent=round(log10(value));
            if abs(exponent)<=3,
                prefixStr=prefixes10(exponent+3);
            else
                prefixStr=prefixes1000(exponent+25);
            end
        end
        function [unit, scalingFactor]=parseUnit(unitString)
            unitString=trim(unitString);
            
        end
    end
end
