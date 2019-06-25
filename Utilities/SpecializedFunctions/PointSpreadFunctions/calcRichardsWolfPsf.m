% [psf psfField]=calcRichardsWolfPsf(xRange,yRange,zRange,wavelength,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture,refractiveIndexOfSample)
%
%       xRange=[-10:.1:10]*1e-6;yRange=[-10:.1:10]*1e-6;zRange=0;
%       objectiveNumericalAperture=0.95;
%       [psf psfField]=calcRichardsWolfPsf(xRange,yRange,zRange,500e-9,1,0,objectiveNumericalAperture,1.0);
% 
function [psf psfField]=calcRichardsWolfPsf(xRange,yRange,zRange,wavelength,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture,refractiveIndexOfSample)
    if nargin<1,
        xRange=0.1*[-10:.2:10]*1e-6;
        yRange=xRange;
        zRange=0;
        wavelength=532e-9;
        pupilFunctorH=1;
        pupilFunctorV=@(X,Y) X*0+0;
        objectiveNumericalAperture=0.95;
    end
    if nargin<8,
        refractiveIndexOfSample=1;
    end
    
    airNA=objectiveNumericalAperture/refractiveIndexOfSample;
    if numel(xRange)>1
        dx=diff(xRange(1:2));
    else
        dx=1;
    end
    if numel(yRange)>1
        dy=diff(yRange(1:2));
    else
        dy=dx;
    end
    k=2*pi/(wavelength/refractiveIndexOfSample);
    
    if ~isa(pupilFunctorH,'function_handle'),
        fl0=1*pupilFunctorH./sqrt(pi*(airNA).^2./(dx*dy));
        calcFieldAtSinglePosition=@(x,y,z) calcDLFieldAtSinglePosition(k,airNA,fl0,x,y,z);
    else
        % f times l_0 as in paper
        fl0=@(sx,sy) pupilFunctorH(sx/airNA,sy/airNA)./sqrt(pi*(airNA).^2./(dx*dy));
        calcFieldAtSinglePosition=@(x,y,z) calcGenericFieldAtSinglePosition(k,airNA,fl0,x,y,z);
    end

    % Iterate over all positions in image space
    psfField=zeros(numel(xRange),numel(yRange),numel(zRange),3);
    for zIdx=1:numel(zRange),
        z=zRange(zIdx);
        for yIdx=1:numel(yRange),
            y=yRange(yIdx);
            for xIdx=1:numel(xRange),
                x=xRange(xIdx);
                psfField(xIdx,yIdx,zIdx,:)=calcFieldAtSinglePosition(x,y,z);
            end
        end
    end
    
    psf=sum(abs(psfField).^2,4);
%     sum(psf(:))
end


function E=calcGenericFieldAtSinglePosition(k,airNA,fl0,x,y,z)
    function result=A(sx,sy,dimIdx)
        ST=sqrt(sx.^2+sy.^2);
        CT=sqrt(1-ST.^2);
        CP=sx./(ST+(ST==0));
        SP=sy./(ST+(ST==0));
        result=repmat(fl0(sx,sy).*sqrt(CT),[1 1 3]).*cat(3,...
            CT+SP.^2.*(1-CT),...
            (CT-1).*CP.*SP,...
            -ST.*CP);
        result=result(:,:,dimIdx);
    end
    integral=zeros(1,1,3);
    for dimIdx=1:3,
        integrand=@(sx,sy,sz) A(sx,sy,dimIdx).*(1./sz).*exp(1i*k*(sx*x+sy*y+sz*z));
        [integral(dimIdx), errorBnd]=quad2d(@(sx,sy) integrand(sx,sy,sqrt(1-(sx.^2+sy.^2))),...
            -airNA,airNA,...
            @(x) -airNA*sqrt(1-(x/airNA).^2),@(x) airNA*sqrt(1-(x/airNA).^2),...
            'AbsTol',1e-6,'MaxFunEvals',5000);
    end
    E=-1i*k/(2*pi)*integral;
end

function E=calcDLFieldAtSinglePosition(k,airNA,fl0,x,y,z)
    rP=sqrt(x.^2+y.^2+z.^2);
    tP=atan2(sqrt(x.^2+y.^2),z);
    phiP=atan2(y,x);
    integrand{1}=@(t) sqrt(cos(t)).*sin(t).*(1+cos(t)).*besselj(0,k*rP.*sin(t).*sin(tP)).*exp(1i*k*rP.*cos(t).*cos(tP));
    integrand{2}=@(t) sqrt(cos(t)).*(sin(t).^2).*besselj(1,k*rP.*sin(t).*sin(tP)).*exp(1i*k*rP.*cos(t).*cos(tP));
    integrand{3}=@(t) sqrt(cos(t)).*sin(t).*(1-cos(t)).*besselj(2,k*rP.*sin(t).*sin(tP)).*exp(1i*k*rP.*cos(t).*cos(tP));
    I=zeros(1,3);
    for idx=1:3,
        I(idx)=quadgk(integrand{idx},0,asin(airNA));
    end
    
    A=k*fl0/2;
    E=zeros(1,1,3);
    E(1)=-1i*A*(I(1)+I(3).*cos(2*phiP));
    E(2)=-1i*A*I(3).*sin(2*phiP);
    E(3)=-2*A*I(2).*cos(phiP);
end