% [pupilFunctionWithoutTiltAndPiston tilt piston]=removeTiltAndPiston(pupilFunction)
%
% Calculates and removes the tilt and piston from a pupil function measurement.
%
% Output parameters:
%     pupilFunctionWithoutTilt: the input with the tilt removed
%     tilt: the removed tilt in units of wavelengths per pixel.
%     piston: the removed piston in units of wavelength.
%
function [pupilFunctionWithoutTiltAndPiston tilt piston]=removeTiltAndPiston(pupilFunction)
    if (nargin<1)
        % Test case
        [X,Y]=ndgrid([-50:50],[-100:100]);
        R=sqrt(X.^2+Y.^2);
        pupilFunction=exp(2i*pi*(0.05*X-0.10*Y  + .001*R.^2 ))*1i.*(R<0.5*max(X(:)));
    end
    inputSize=size(pupilFunction);
    
    [output pupilFunctionWithoutTiltAndPiston]=dftregistration(ones(inputSize),ifftshift(pupilFunction),100);
    piston=-output(2)/(2*pi);
    tilt=output(3:4)./inputSize;
    
    pupilFunctionWithoutTiltAndPiston=fftshift(pupilFunctionWithoutTiltAndPiston);
end