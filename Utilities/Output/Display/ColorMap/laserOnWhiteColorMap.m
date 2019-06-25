% RGB=laserOnWhiteColorMap(nbEntries,wavelength,powerLevel)
%
% Creates a color map to be used with the command colormap
% All arguments are optional.
%
% Example usage:
%    figure;
%    image([0:.001:1]*255);
%    colormap(laserOnWhiteColorMap(256,532e-9,1))
function RGB=laserOnWhiteColorMap(nbEntries,wavelength,powerLevel)
    if (nargin<1 || isempty(nbEntries))
        nbEntries=64;
    end
    if (nargin<2 || isempty(wavelength))
        wavelength=633e-9;
    end
    if (nargin<3)
        powerLevel=1;
    end

     intensities=([1:nbEntries]-1)/max(1,nbEntries-1);
     RGB=repmat(1-intensities.',[1 3])+repmat(spectrumToRGB(wavelength,powerLevel/3,false),[nbEntries 1]);
     
     %Saturate the 'eye'
     RGB=RGB-repmat(min(0,min(RGB,[],2)),[1 3]);
     RGB=min(1,RGB);
     
%      RGB(1,:)=1; % Make sure that zero intensity is white
end