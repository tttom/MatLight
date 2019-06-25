% Used for light sheet simulations and deconvolution
function psf=calcDetectionPsf(xRange,yRange,zRange,u,v,detection,refractiveIndexOfSample)    
    minYSize=512;
    if (length(yRange)>=minYSize)
        extendedYRange=yRange;
    else
        %Crank up the resolution temporarily to avoid under sampling of the pupil
        if (length(yRange)>1)
            diffY=diff(yRange(1:2));
        else
            diffY=diff(xRange(1:2));
        end
        extendedYRange=diffY*([1:minYSize]-(floor(minYSize/2)+1))+yRange(floor(end/2)+1);
    end
    
    pupilFunctor=@(U,V) exp(2i*pi*(u*U+v*V));
%  old code, lots of approximations:   psf=calcPsf(xRange,extendedYRange,zRange,detection,pupilFunctor,refractiveIndexOfSample); 
    numericalApertureInSample=detection.objective.numericalAperture/detection.objective.refractiveIndex; %Because this is an air objective, so the NA is the same in the sample
    %Assume circular polarization propagating along the x-axis
    psf=calcVectorialPsf(xRange,yRange,zRange,detection.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalApertureInSample,refractiveIndexOfSample,detection.magnification,detection.tubeLength);
        
    if (length(extendedYRange)~=length(yRange))
        psf=psf(:,floor(end/2)+1+[1:length(yRange)]-(floor(length(yRange)/2)+1),:);
    end
end