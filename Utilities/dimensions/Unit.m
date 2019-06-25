% 
% A class bundling all kinds of units
%
classdef Unit
    properties (Constant)
        % metric prefixes first
        y=1e-24; yocto=Unit.y;
        z=1e-21; zepto=Unit.z;
        a=1e-18; atto=Unit.a;
        f=1e-15; femto=Unit.f;
        p=1e-12; pico=Unit.p;
        n=1e-9; nano=Unit.n;
        u=1e-6; micro=Unit.u;
        m=1e-3; milli=Unit.m;
        
        c=1e-2; centi=Unit.c;
        d=1e-1; deci=Unit.d;
        da=1e1; deca=Unit.da;
        h=1e2; hecto=Unit.h;
        
        k=1e3; kilo=Unit.k;
        M=1e6; mega=Unit.M;
        G=1e9; giga=Unit.G;
        T=1e12; tera=Unit.T; 
        P=1e15; peta=Unit.P;
        E=1e18; exa=Unit.E;
        Z=1e21; zetta=Unit.Z;
        Y=1e24; yotta=Unit.Y;
        
        % Specific to volume
        L=0.1^3; liter=Unit.L % m^3
        
        % English units
        imperial=Imperial();
        us=US();
    end
end
