%Adapted from Michael's code
%
% numericalAperture is in air!
function psf=calcBesselPsf(xRange,yRange,zRange,lambda,numericalApertureInAir,refractiveIndexOfSample)
    c=299792458;
    pol=[1 0]; % Polarization
    
    [Y,X,Z]=meshgrid(yRange,xRange,zRange);
    
    [T,R,Z2D]=cart2pol(X(:,:,1),Y(:,:,1),Z(:,:,1)); %Only the phase of the field changes in the propagation direction.
    R=max(R,eps(1));
    omega=2*pi*c/(lambda/refractiveIndexOfSample);
    gamma=asin(numericalApertureInAir/refractiveIndexOfSample);
    lNumber=0; % l-number
    tilt=0;
    
    erho=	(1i/2.*exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*(2.*pol(1).*c.*lNumber.*besselj(lNumber,(R.*omega.*sin(gamma))./c) + pol(2).*R.*omega.*(besselj(-1 + lNumber,(R.*omega.*sin(gamma))./c) - besselj(1 + lNumber,(R.*omega.*sin(gamma))./c)).*cos(gamma).*sin(gamma)))./(c.*R);	
    ephi=	-(exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*(2.*pol(2).*c.*lNumber.*besselj(lNumber,(R.*omega.*sin(gamma))./c).*cos(gamma) + pol(1).*R.*omega.*(besselj(-1 + lNumber,(R.*omega.*sin(gamma))./c) - besselj(1 + lNumber,(R.*omega.*sin(gamma))./c)).*sin(gamma)))./(2*c.*R);	
    ez=	(pol(2).*exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*omega.*besselj(lNumber,(R.*omega.*sin(gamma))./c).*sin(gamma).^2)./c;	
    hrho=	(exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*(2.*c.*lNumber.*besselj(lNumber,(R.*omega.*sin(gamma))./c).*(pol(2) - pol(1).*cos(gamma)) + pol(1).*R.*omega.*besselj(-1 + lNumber,(R.*omega.*sin(gamma))./c).*sin(2*gamma)))./(2*c.*R);	
    hphi=	((-1i).*exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*(c.*lNumber.*besselj(lNumber,(R.*omega.*sin(gamma))./c).*(pol(2) - pol(1).*cos(gamma)) - pol(2).*R.*omega.*besselj(-1 + lNumber,(R.*omega.*sin(gamma))./c).*sin(gamma)))./(c.*R);	
    hz=	((-1i).*pol(1).*exp((1i.*(c.*lNumber.*T - c.*tilt.*omega + Z2D.*omega.*cos(gamma)))./c).*omega.*besselj(lNumber,(R.*omega.*sin(gamma))./c).*sin(gamma).^2)./c;	

    efx=erho.*cos(T)-ephi.*sin(T);
    efy=erho.*sin(T)+ephi.*cos(T);
    efz=ez;
    hfx=hrho.*cos(T)-hphi.*sin(T);
    hfy=hrho.*sin(T)+hphi.*cos(T);
    hfz=hz;
    %%%%%
    psf=abs(efx).^2+abs(efy).^2+abs(efz).^2;
    psf=repmat(psf,[1 1 size(Z,3)]); %Extend the beam along the optical axis.
    
    psf=psf./sum(sum(psf(:,:,1))); %Normalize cross section intensity
end