% psf=calcPsf(xRange,yRange,zRange,objective,functor,refractiveIndexOfSample,projectionDimension)
%
% Calculates the 3D intensity point spread function at the grid specified by
% x/y/zRange for a pupil function given by functor(U,V) where U and V are
% normalized carthesian pupil coordinates, in a medium with
% refractive index refractiveIndexOfSample and an objective specified by
% the struct objecive with fields:
%    -wavelength
%    -numericalApertureInAir
%    -magnification
%    -tubeLength
%
function psf=calcPsf(xRange,yRange,zRange,objective,functor,refractiveIndexOfSample,projectionDimension)
    if (nargin<7)
        projectionDimension=[];
    end
    
    gridSize=[length(xRange) length(yRange) length(zRange)];
    
    targetGridSize=gridSize;
    for (projDimIdx=projectionDimension)
        targetGridSize(projDimIdx)=1;
    end
    
    partSize=ceil(64e6/prod(gridSize(1:2)));
    
    focalLengthInSample=(objective.tubeLength/objective.magnification)*refractiveIndexOfSample;
    objectiveSinMaxAngleInSample=objective.numericalApertureInAir/refractiveIndexOfSample;
    apertureRadiusInSample=focalLengthInSample*objectiveSinMaxAngleInSample;
    wavelengthInSample=objective.wavelength/refractiveIndexOfSample;
    
    cutOffSpatialFrequency=2*objective.numericalApertureInAir/objective.wavelength;
    
    if (length(xRange)>1)
        dx=diff(xRange(1:2));
    else
        dx=1;
    end
    if (length(yRange)>1)
        dy=diff(yRange(1:2));
    else
        dy=1;
    end
    decenterInWaves=[xRange(floor(end/2)+1) yRange(floor(end/2)+1)]*cutOffSpatialFrequency/2;
    uRange=2*([1:gridSize(1)]-1-floor(gridSize(1)/2))/(gridSize(1)*dx)/cutOffSpatialFrequency;
    vRange=2*([1:gridSize(2)]-1-floor(gridSize(2)/2))/(gridSize(2)*dy)/cutOffSpatialFrequency;
    pupilCenterToFocus=sqrt(focalLengthInSample^2-apertureRadiusInSample^2);
    wRange=(sqrt((pupilCenterToFocus+zRange).^2+apertureRadiusInSample^2)-(focalLengthInSample+zRange))/wavelengthInSample;
    
    [U,V]=ndgrid(uRange,vRange);
    R2=U.^2+V.^2;
    
    circularPupilTransmission=ones(1,class(uRange))*(R2<1);
    circularPupilTransmission=circularPupilTransmission./sqrt(sum(circularPupilTransmission(:)));
    pupil2D=functor(U,V).*circularPupilTransmission.*exp(2i*pi*(decenterInWaves(1)*U+decenterInWaves(2)*V));
    clear U V;
    
    psf=zeros(targetGridSize,class(xRange));
    %Do block by block transforms to make optimal use of memory and CPU
    for (partIdx=1:partSize:numel(wRange))
        thisPartSize=min(partSize,length(wRange)-partIdx+1);
        pupil=R2(:)*wRange(partIdx-1+[1:thisPartSize]);
        pupil=reshape(exp(2i*pi*pupil),[gridSize(1:2) thisPartSize]);
        pupil=repmat(pupil2D,[1 1 thisPartSize]).*pupil;

        psfPart=fftshift(fftshift(ifft2(pupil),1),2)*sqrt(numel(pupil2D)); %Center ifftshift not required because of abs value in the next step
        psfPart=abs(psfPart).^2;

        for (projDimIdx=projectionDimension)
            psfPart=sum(psfPart,projDimIdx);
        end
    
        psf(1:size(psfPart,1),1:size(psfPart,2),partIdx-1+[1:size(psfPart,3)])=psfPart;
    end
end