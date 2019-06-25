% [psf, psfField, varargout] = calcVectorialPsf(xRange,yRange,zRange,wavelength,...
%                                pupilFunctorH,pupilFunctorV,...
%                                objectiveNumericalAperture,refractiveIndexOfSample,...
%                                projectionDimensions)
%
% Calculates the 3D point spread function at the grid specified by
% x/y/zRange for a pupil function given by pupilFunctorH/V(U,V) where U and V
% are normalized cartesian pupil coordinates, in a medium with
% refractive index refractiveIndexOfSample and an objective. The horizontal
% polarization given by pupilFunctorH is along the first dimension, and the
% vertical given by pupilFunctorV is along the second dimension of the
% output.
%
% x/y/zRange are the range of the image produced
%
% This function gives approximately the same results as PSFLab:
% http://onemolecule.chem.uwm.edu/software , though it is significantly
% faster. Please contact Tom Vettenburg in case you suspect discrepancies
% with the theoretical or simply for more information.
%
% Output:
%     psf: the single photon intensity
%     psfField: the electric field components (x,y,z listed in the fourth
%               dimension of the matrix). For non-vectorial calculations
%               (pupilFunctorV==[]), this matrix is of a maximum of three
%               dimensions.
%     varargout: the higher order nonlinear intensities. This is
%     effectivelly the same as psf.^(N-1) unless projectionDimensions is
%     specified.
%
% Inputs:
%     x/y/zRange: The sample position in metric cathesian coordinates.
%                 Singleton dimensions x and y are treated for normalization
%                 as if they were Nyquist sampled.
%     wavelength: The wavelength in the same metric units.
%     pupilFunctorH: A function returning the complex horizontal pupil field (along X) 
%                    as a function of carthesian coordinates normalized to the pupil radius.
%                    When a scalar is specified, a constant field of this value is assumed for the whole circular pupil.
%                    Default: unity transmission inside the pupil
%     pupilFunctorV: A function returning the complex vertical pupil field (along Y)
%                    as a function of carthesian coordinates normalized to the pupil radius.
%                    When a scalar is specified, a constant field of this value is assumed for the whole circular pupil.
%                    When a vector or a matrix is specified, the values
%                    will be tiled uniformely over the 2x2 square covering
%                    the pupil disk.
%                    A scalar calculation will be done instead when nothing
%                    or the empty list is given. Specify 0 when only
%                    horizontal input fields are required. Don't use the
%                    empty matrix unless scalar calculations are required.
%                    Default: []: scalar calculation.
%     numericalApertureInSample: (default 1.0) The objective's numerical aperture (including refractive index, n, of the medium: n*sin(theta)).
%     refractiveIndexOfSample: (default 1.0 for vacuum) The refractive index of the medium at the focus.
%                    Cover slip correction is assumed in the calculation, hence
%                    this only scales the sample grid.
%     projectionDimensions:
%               -when omitted or empty ([]), the full data cube is returned
%               -when a vector with integers, the data is integrated along
%                        the dimensions indicated in the vector.
%                   Default: no projection ([])
%
% Examples:
%     % Plot a 1D cross section through the PSF of circularly polarized beam send of wavelength 532nm, focussed by an objective with NA=0.42 in water. 
%     xRange=[-10:.1:10]*1e-6;
%     psf=calcVectorialPsf(xRange,0,0,532e-9,@(U,V) 1,@(U,V) 1i,0.42,1.33);
%     figure(); plot(xRange*1e6,psf);
%
%     xRange=[-10:.1:10]*1e-6;yRange=[-10:.1:10]*1e-6;zRange=[-10:.1:10]*1e-6;
%     objectiveNumericalAperture=asin(1./sqrt(2));
%     pupilFunctor=@(U,V) sqrt(U.^2+V.^2)>.9; %Bessel beam with 10% open fraction
%     [psf psfField]=calcVectorialPsf(xRange,yRange,zRange,500e-9,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,1.0);
%     psfProj=calcVectorialPsf(xRange,yRange,zRange,500e-9,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,1.0,[2]);
%
%     img=squeeze(psf(:,1+floor(end/2),:)).';
%     img=img./repmat(mean(img,2),[1 size(img,2)]);
%     figure();
%     imagesc(xRange*1e6,zRange*1e6,img);axis equal;xlabel('x [\mu m]');ylabel('z [\mu m]');
%
function [psf, psfField, varargout]=calcVectorialPsf(xRange,yRange,zRange,wavelength,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture,refractiveIndexOfSample,projectionDimensions)
    debug = false;
    if (nargin<1 || isempty(xRange))
        xRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<2 || isempty(yRange))
        yRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<3 || isempty(zRange))
        zRange=[-5000:500:5000]*1e-9;
    end
    if (nargin<4 || isempty(wavelength))
        wavelength=532e-9;
    end
    if (nargin<5 || isempty(pupilFunctorH))
        % Don't change the following, it is a sensible default in case the user wants Vertical polarization only!
        pupilFunctorH=@(normalU,normalV) 0.0; % Polarization along the first dimension
    end
    if (nargin<6)
        pupilFunctorV=[]; % Polarization along the second dimension, default: only along the first dimension
    end
    if (nargin<7 || isempty(objectiveNumericalAperture))
        objectiveNumericalAperture=1.0;
    end
    if (nargin<8 || isempty(refractiveIndexOfSample))
        refractiveIndexOfSample=1.0;
    end
    if (nargin<9)
        projectionDimensions=[];
    end
    
    calculationSize=[numel(xRange) numel(yRange) numel(zRange)]; % the size without projection, bar field coordinates
    outputSize=calculationSize;
    outputSize(projectionDimensions)=0; % the size of the spatial coordinates of the output data set
    
    % Check how many multi-photon orders of the intensity have to be calculated
    highestOrderIntensityRequired = max(1, nargout-1);
    
    vectorialCalculation=~isempty(pupilFunctorV);
    if (debug && ~vectorialCalculation)
        logMessage('Starting a scalar calculation of the PSF...');
    end
    
    % convert the horizontal pupil to a function if needed
    if ischar(pupilFunctorH)
        % convert specific polarization states to functions
        pupilFunctorV=pupilFunctorH; % let the vertical component handle it
        pupilFunctorH=@(U,V) 1;
    else
        pupilFunctorH=functionalizePupilDescription(pupilFunctorH);
    end
    % convert the vertical pupil to a function if needed
    if ischar(pupilFunctorV)
        % convert specific polarization states to functions
        switch lower(pupilFunctorV)
            case {'horizontal','hori'}
                pupilFunctorV=0;
            case {'vertical','vert'}
                pupilFunctorV = pupilFunctorH;
                pupilFunctorH = 0;
            case 'left'
                pupilFunctorV = @(U,V) 1i*pupilFunctorH(U,V);
            case 'right'
                pupilFunctorV = @(U,V) -1i*pupilFunctorH(U,V);
        end
    else
        % Convert numeric input to interpolating functions
        pupilFunctorV = functionalizePupilDescription(pupilFunctorV);
    end
        
    objectiveSinMaxHalfAngleInSample=objectiveNumericalAperture/refractiveIndexOfSample;
    
    % Determine the requested step size for each non-singleton dimension
    if all(calculationSize>0)
        sampleDelta(1) = abs(diff(xRange([1 min(2,end)])));
        sampleDelta(2) = abs(diff(yRange([1 min(2,end)])));
        sampleDelta(3) = abs(diff(zRange([1 min(2,end)])));
    
        % Setting constants
        minPupilSize = [1 1]*256; % To ensure that rapid pupil changes are properly sampled
        maxPupilSize = [1 1]*1024; % To avoid memory problems

        %The minimum pupil size to avoid PSF replication
        requiredPupilSize = 2*ceil([numel(xRange) numel(yRange)].*sampleDelta(1:2)/(wavelength/(objectiveSinMaxHalfAngleInSample*refractiveIndexOfSample)));
        if (debug &&  any(requiredPupilSize>maxPupilSize))
            logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle.\nThis will cause replicas.',[maxPupilSize,requiredPupilSize]);
        end

        %Check if pupil sampling not too sparse for the defocus we intend to simulate
        maxSampleNA = objectiveSinMaxHalfAngleInSample*(1-0.25/(max(maxPupilSize)/2)); % Use of 'max' because the sampling rate near the edge for NA=1 diverges
        minPupilSizeToHandleDefocus = minPupilSize+[1 1]*max(abs(zRange/(wavelength/refractiveIndexOfSample)))*4*objectiveSinMaxHalfAngleInSample*maxSampleNA/sqrt(1-maxSampleNA^2);
        if any(minPupilSize<minPupilSizeToHandleDefocus)
            if debug
                logMessage('A minimum pupil size of (%0.0f,%0.0f) is required to handle the specified defocus.',minPupilSizeToHandleDefocus);
            end
            requiredPupilSize=max(requiredPupilSize,minPupilSizeToHandleDefocus);
        end
        if (debug && any(minPupilSizeToHandleDefocus>maxPupilSize))
            logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle. This can cause aliasing artefacts.',[maxPupilSize,minPupilSizeToHandleDefocus]);
        end

        pupilGridSize = ceil(min(maxPupilSize, max(minPupilSize,requiredPupilSize)));
        if debug
            logMessage('Setting the pupil sampling grid size to (%0.0f,%0.0f)',pupilGridSize);
        end

        %Check if CZT of FFT is faster
        nextFastestArray = [2 2 3 5 5 7 7 9 9 10 12 12 14 14 16 16 20 20 20 20 22 22 25 25 25 28 28 28 32 32 32 32 33 36 36 36 39 39 39 40 42 42 48 48 48 48 48 48 52 52 52 52 54 54 64 64 64 64 64 64 64 64 64 64 66 66 72 72 72 72 72 72 80 80 80 80 80 80 80 80 84 84 84 84 96 96 96 96 96 96 96 96 96 96 96 96 98 98 100 100 110 110 110 110 110 110 110 110 110 110 120 120 120 120 120 120 120 120 120 120 126 126 126 126 126 126 128 128 130 130 132 132 140 140 140 140 140 140 140 140 144 144 144 144 150 150 150 150 150 150 160 160 160 160 160 160 160 160 160 160 162 162 176 176 176 176 176 176 176 176 176 176 176 176 176 176 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 196 196 196 196 200 200 200 200 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 224 240 240 240 240 240 240 240 240 240 240 240 240 240 240 240 240 256 256 256 256 256 256 256 256 256 256 256 256 256 256 256 256 260 260 260 260 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 320 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 352 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 384 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 416 416 416 416 416 416 416 416 416 416 416 416 416 416 416 416 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 448 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 512 520 520 520 520 520 520 520 520 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 560 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 640 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 768 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 800 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 832 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 896 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1024 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1040 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1120 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1280 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1296 1300 1300 1300 1300 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1320 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1344 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1375 1386 1386 1386 1386 1386 1386 1386 1386 1386 1386 1386 1400 1400 1400 1400 1400 1400 1400 1400 1400 1400 1400 1400 1400 1400 1404 1404 1404 1404 1408 1408 1408 1408 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1452 1456 1456 1456 1456 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1584 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1600 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1620 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1664 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1680 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1728 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 1792 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2048 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2080 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2240 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2520 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2560 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2592 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2640 2646 2646 2646 2646 2646 2646 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2688 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2816 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 2912 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3080 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3120 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3328 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3380 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3520 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3584 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3640 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 3840 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096];
        nextFastest = @(x) nextFastestArray(min(end,ceil(x)));

        cztFftSize = nextFastest(pupilGridSize + calculationSize(1:2)-1);
        Dt = (wavelength/(2*objectiveNumericalAperture))./sampleDelta(1:2); % The support size required in the pupil plane
        fftFftSize = nextFastest(pupilGridSize.*Dt);
        workCZT = 2*(cztFftSize.*log2(cztFftSize));
        workFFT = fftFftSize.*log2(fftFftSize);
        if debug
            fftPupilSamples = fftFftSize./Dt;
            cztRelTimes = workCZT./workFFT;
            cztPupilSamples = cztFftSize-(calculationSize(1:2)-1);
            cztRelAccuracy = cztPupilSamples./fftPupilSamples;
            cztRelToCurrentAccuracy = cztPupilSamples./pupilGridSize;
            fftRelToCurrentAccuracy = fftPupilSamples./pupilGridSize;
            if prod(cztRelTimes)<1
                logMessage('CZT method will save %0.0f%% (%0.0f%%x%0.0f%%) of time over FFT.',100*(1-[prod(cztRelTimes) cztRelTimes]));
            else
                logMessage('FFT method will save %0.0f%% (%0.0f%%x%0.0f%%) of time.',100*(1-1./[prod(cztRelTimes) cztRelTimes]));
            end
            if prod(cztRelAccuracy)>1
                logMessage('CZT could increase sampling density by %0.0f%% (%0.0f%%x%0.0f%%) over FFT, and by %0.0f%% (%0.0f%%x%0.0f%%) over current.',...
                    100*([prod(cztRelAccuracy) cztRelAccuracy prod(cztRelToCurrentAccuracy) cztRelToCurrentAccuracy]-1));
            else
                logMessage('FFT could increase sampling density by %0.0f%% (%0.0f%%x%0.0f%%) over CZT, and by %0.0f%% (%0.0f%%x%0.0f%%) over current.',...
                    100*(1./[prod(cztRelAccuracy) cztRelAccuracy prod(fftRelToCurrentAccuracy) fftRelToCurrentAccuracy]-1));
            end
        end
        % Adapt pupil sampling if this can gain speed, cztw() will pick the fastest algorithm
        pupilDimensions = 2.0*[1 1]; % The area covered by the grid of size pupilGridSize in normalized pupil coordinates.
        if any(workFFT<workCZT) && false % FFT processing disabled until it is fixed in cztw
            % FFT will be faster than CZT
            if debug
                logMessage('Increasing sampling density to match FFT.');
            end
            % Change the sampling size to match that for an FFT
            pupilGridSize(workFFT<workCZT) = fftFftSize(workFFT<workCZT)./Dt(workFFT<workCZT);

            % Only calculate points that are inside the unit radius aperture using floor.
            % Change the dimensions of the pupil accordingly:
            pupilDimensions = pupilDimensions.*floor(pupilGridSize)./pupilGridSize;
            pupilGridSize = floor(pupilGridSize);
        end

        %Choose the pupil grid
        wavelengthInSample = wavelength/refractiveIndexOfSample;
        uRange = objectiveSinMaxHalfAngleInSample*2*[-floor(pupilGridSize(1)/2):floor((pupilGridSize(1)-1)/2)]/pupilGridSize(1);
        vRange = objectiveSinMaxHalfAngleInSample*2*[-floor(pupilGridSize(2)/2):floor((pupilGridSize(2)-1)/2)]/pupilGridSize(2);

        [U,V] = ndgrid(uRange,vRange);
        sinApAngle2 = U.^2+V.^2;
        apertureFieldTransmission = double(sinApAngle2<objectiveSinMaxHalfAngleInSample^2);
        totalPowerForUnityIrradiation = numel(U)*2*pi*(1-sqrt(1-objectiveSinMaxHalfAngleInSample^2))/(4*objectiveSinMaxHalfAngleInSample^2); % integral of unit pupil squared (accounting for the sine condition correction) / total 2D integral = 2pi(1-sqrt(1-NA^2))/ 4 NA^2
        sinApAngle = sqrt(apertureFieldTransmission.*sinApAngle2);
        cosApAngle = apertureFieldTransmission.*sqrt(1-apertureFieldTransmission.*sinApAngle2);
        clear sinApAngle2;
        %Scale so that the total intensity is 1 for a unity uniform illumination
        apertureFieldTransmission = apertureFieldTransmission./sqrt(totalPowerForUnityIrradiation); % sqrt because the field is the sqrt of the intensity which is preserved

        % Instead of the Herschel condition, use the sine condition.
        % The field on the back-aperture hemisphere should be scaled with sqrt(cos(t)) to account for its curvature.
        % The change from spherical to Cartesian integration variables: 1/cos(t)
        sineConditionFieldScaling = (cosApAngle+(cosApAngle==0)).^-0.5; % avoid division by zero by setting 1/sqrt(0) to 1/sqrt(1)
        apertureFieldTransmission = sineConditionFieldScaling.*apertureFieldTransmission;

        if debug
            logMessage('Total intensity in = %0.3f',sum(abs(apertureFieldTransmission(:)).^2));
        end

        pupilFunctionX = apertureFieldTransmission.*pupilFunctorH(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
        if vectorialCalculation
            pupilFunctionY = apertureFieldTransmission.*pupilFunctorV(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
            %Convert pupil function to polar coordinates
            T = atan2(V,U); CT = cos(T); ST = sin(T);
            pupilFunctionR =  CT.*pupilFunctionX + ST.*pupilFunctionY; % Radial component is rotated by the focusing
            pupilFunctionA = -ST.*pupilFunctionX + CT.*pupilFunctionY; % Azimutal component is unaffected by the focusing
            %Calculate the polarization change due to focussing
            pupilFunctionZ = sinApAngle.*pupilFunctionR;
            pupilFunctionR = cosApAngle.*pupilFunctionR;
            %Convert back to Cartesian coordinates
            pupilFunctionX = CT.*pupilFunctionR - ST.*pupilFunctionA;
            pupilFunctionY = ST.*pupilFunctionR + CT.*pupilFunctionA;
            clear pupilFunctionR pupilFunctionA CT ST T;

            pupilFunction2D = cat(3,pupilFunctionX,pupilFunctionY,pupilFunctionZ);
            clear pupilFunctionX pupilFunctionY pupilFunctionZ;
        else
            pupilFunction2D = pupilFunctionX;
            clear pupilFunctionX;
        end
        clear apertureFieldTransmission sinApAngle cosApAngle sineConditionFieldScaling;

        %Calculate the focal plain fields
        [psfField, psfIntensities] = czt2andDefocus(pupilFunction2D,pupilDimensions,objectiveSinMaxHalfAngleInSample,xRange/wavelengthInSample,yRange/wavelengthInSample,zRange/wavelengthInSample, projectionDimensions, highestOrderIntensityRequired);
        
    else
        psfField = zeros([outputSize 3],class(xRange));
        psfIntensities = zeros([outputSize highestOrderIntensityRequired],class(xRange));
    end

    % Rename the output
    psf = psfIntensities(:,:,:,1);
    if size(psfIntensities,4) > 2
        varargout = mat2cell(psfIntensities(:,:,:,2:end),size(psfIntensities,1),size(psfIntensities,2),size(psfIntensities,3),ones(1,size(psfIntensities,4)-1));
    else
        if size(psfIntensities,4) == 2
            varargout = {psfIntensities(:,:,:,2)};
        end
    end
%     clear psfIntensities;
end

% Calculate the partial spectrum of x using the chirp z transform.
% This returns the complex field.
%
% x is the pupil and should not be ifftshifted, implicit zero padding to the right!
% pupilDimensions: should be [2 2] in case of CZT
% objectiveSinMaxHalfAngleInSample: the input matrix must cover this disk exactly.
% the following arguments, also specifiable as a list, are:
%     nxRange and nyRange: the sample points in the units of wavelength in the sample medium.
%     nzRange: the sample points in the z dimension in units of wavelength in the sample.
%     projectionDimensions: (optional, default none) The dimension along which an integration is done
%     highestOrderIntensityRequired: (optional, default depends on nargout) If specified, (higher order) intensities upto this number are returned as well.
%
% If more than one output argument is specified, the first and higher order
% intensities will be returned as 3D arrays stacked into a single 4D array.
%
function [f, psfIntensities] = czt2andDefocus(x,pupilDimensions,objectiveSinMaxHalfAngleInSample,nxRange,nyRange,nzRange, projectionDimensions, highestOrderIntensityRequired)
    if (nargin<6)
        projectionDimensions = [];
    end
    if (nargin<7)
        highestOrderIntensityRequired = max(0,nargout-1);
    end
    
    normalizedMaxFreqInSample = objectiveSinMaxHalfAngleInSample*pupilDimensions;
    
    %Prepare the output matrix with zeros
    inputSize = size(x);
    %inputSize = x(1)*0+inputSize;%Cast to same class as the x input
    outputSize = [length(nxRange) length(nyRange) length(nzRange), inputSize(3:end)]; % Add dimension
    if (~isempty(projectionDimensions))
        outputSize(projectionDimensions) = 1;
    end
    f = zeros(outputSize,class(x));
    psfIntensities = zeros([outputSize(1:3) highestOrderIntensityRequired],class(x));
    
    uRange = normalizedMaxFreqInSample(1)*[-floor(inputSize(1)/2):floor((inputSize(1)-1)/2)]/inputSize(1);
    vRange = normalizedMaxFreqInSample(2)*[-floor(inputSize(2)/2):floor((inputSize(2)-1)/2)]/inputSize(2);
    [U,V] = ndgrid(uRange,vRange);
    clear uRange vRange;
    R2 = min(1.0,U.^2+V.^2);
    clear U V;
    cosHalfAngleInSampleMatrix = sqrt(1-R2);
    clear R2;
    
    %Loop through the z-stack
    pupil = zeros(size(x));
    for zIdx = 1:numel(nzRange)
        normalizedZ = nzRange(zIdx);
        phaseOffsetOfAxialWaveletInRad = 2*pi*mod(normalizedZ,1); %Center on focal point
        % Calculate the phase delay due to z-displacement with respect to
        % the wavelet leaving the center of the pupil. This causes the Gouy
        % phase shift.
        unityDefocusInRad = 2*pi*(cosHalfAngleInSampleMatrix-1); %Assume that virtual image of aperture at infinity
        relativeDefocusPhaseDelayAcrossPupilInRad = normalizedZ*unityDefocusInRad;

        pupil(1:size(x,1),1:size(x,2),:) = -1i*x.*repmat(exp(1i*(phaseOffsetOfAxialWaveletInRad+relativeDefocusPhaseDelayAcrossPupilInRad)),[1 1 inputSize(3:end)]);
        % The -1i is just to get the absolute field phase correct in the output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psfSlice = cztFromRanges(pupil,nxRange*normalizedMaxFreqInSample(1),nyRange*normalizedMaxFreqInSample(2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Project the output before continuing to save memory
        psfSliceIntensity = sum(abs(psfSlice).^2,3);
        if highestOrderIntensityRequired>1
            psfSliceIntensity = repmat(psfSliceIntensity,[1 1 highestOrderIntensityRequired]);
            for photonNb = 2:highestOrderIntensityRequired
                psfSliceIntensity(:,:,photonNb) = psfSliceIntensity(:,:,photonNb).^photonNb;
            end
        end
        if ~isempty(projectionDimensions)
            for projIdx=1:size(projectionDimensions, 2)
                projectionDimension=projectionDimensions(projIdx);
                if (projectionDimension>=3)
                    projectionDimension=projectionDimension-1;
                end
                psfSlice=sum(psfSlice,projectionDimension);
                psfSliceIntensity=sum(psfSliceIntensity,projectionDimension);
            end
        end
        if any(projectionDimensions == 3)
            f(:,:,1,:) = f(:,:,1,:) + psfSlice;
            psfIntensities(:,:,1,:) = psfIntensities(:,:,1,:) + psfSliceIntensity;
        else
            f(:,:,zIdx,:) = psfSlice;
            psfIntensities(:,:,zIdx,:) = psfSliceIntensity;
        end
    end % of for loop over zIdx
end


function pupilFunctor = functionalizePupilDescription(pupilDescription)
    pupilFunctor=pupilDescription; % return the same if input cannot be handled
    %When scalars are specified instead of function, assume constant input fields
    if isa(pupilDescription,'function_handle')
        pupilFunctor=pupilDescription;
    else if ~isempty(pupilDescription)
            if isscalar(pupilDescription)
                pupilFunctor=@(normalU,normalV) pupilDescription;
            else
                % if a matrix or vector is given, tile the values across
                % the 2x2 square covering the pupil
                if isnumeric(pupilDescription) && ndims(pupilDescription)<3
                    descriptionMatrixSize=size(pupilDescription); descriptionMatrixSize(end+1:2)=1;
                    xRangeOfDescriptionMatrix=2*([1:descriptionMatrixSize(1)]-.5)./descriptionMatrixSize(1)-1;
                    yRangeOfDescriptionMatrix=2*([1:descriptionMatrixSize(2)]-.5)./descriptionMatrixSize(2)-1;
                    if numel(xRangeOfDescriptionMatrix)>1 && numel(yRangeOfDescriptionMatrix)>1
                        pupilFunctor=@(normalU,normalV) interp2(xRangeOfDescriptionMatrix,yRangeOfDescriptionMatrix,pupilDescription,normalU,normalV,'*nearest',0);
                    else
                        if numel(xRangeOfDescriptionMatrix)<=1
                            pupilFunctor=@(normalU,normalV) interp1(yRangeOfDescriptionMatrix,pupilDescription,normalV,'*nearest','extrap');
                        else
                            pupilFunctor=@(normalU,normalV) interp1(xRangeOfDescriptionMatrix.',pupilDescription,normalU,'*nearest','extrap');
                        end
                    end
                end
            end
        end
    end
end