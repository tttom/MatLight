% F = cztFromRanges(T, varargin)
%
% E.g. F = cztFromRanges(T, xfRange, yfRange)
%
% Inputs:
%   varargin: should be linearly increasing vectors.
%
% Output:
%   F: The spectrum at the points indicated by the ranges
%
% Calculate the partial spectrum of T using the 2-dimensional chirp z transform.
% This would be the same as F = fftshift(fft2(ifftshift(T)))*sqrt(numel(T)/diff(xfRange(1:2))) for:
% xfRange = [0:size(T,1)-1]-floor(size(T,1)/2); yfRange = [0:size(T,2)-1]-floor(size(T,2)/2);
% The L2-norm is scaled so that it would equal the L2-norm of the input for large output sizes.
%
% When the input frequencies are different from [0:size(T,1)-1]-floor(size(T,1)/2),
% the outputRanges can be converted as:
% relOutputRange = (outputRange - inputRange(1+floor(end/2)))./diff(inputRange(1:2));
%
% F = fftshift(fft(ifftshift(T)));
% F = cztFromRanges(T, fRange*numel(tRange))*sqrt(numel(functionT));
%     
% The units of x/yfRange are thus twice the highest possible spatial frequency
% of the field PSF when the pupil disc fits the sample area of input T.
% Nyquist sampling is done for x/yfRange = [from:0.5:to];
%
% T should not be ifftshifted, implicit zero padding to the right!
%
function F = cztFromRanges(T, varargin)
    if nargin<1
        tRange = (-500:499).';
        T = sin(2*pi*3*tRange./numel(tRange));
    end
    if nargin<2
        fRange = .1*(-500:499).';
        varargin{1} = fRange;
    end

    inputSize = size(T);
    outputRanges = varargin;
    nbRanges = numel(outputRanges);
    
    deltaRng = zeros(1,nbRanges); % Vector with the step size in each range (default 0, may be updated later)
    M = zeros(1,nbRanges); A = M;
    for (dimIdx = 1:nbRanges)
        rng = outputRanges{dimIdx};
        if ~isempty(rng)
            M(dimIdx) = numel(rng); % The output length of the transform
            A(dimIdx) = exp(2i*pi*rng(1+floor(end/2))/inputSize(dimIdx)); % Offset on contour (fftshifted, hence~/2)
            if (M(dimIdx)>1)
                deltaRng(dimIdx) = diff(rng(1:2));
                % else 0
            end
        else
            % Dimension ignored
            M(dimIdx) = inputSize(dimIdx);
            A(dimIdx) = NaN;
        end
    end
    deltaOmegaInCycles = deltaRng./inputSize(1:nbRanges);
    W = exp(-2i*pi*deltaOmegaInCycles); % The ratio between the points on the spiral contour
    
    if (any(deltaOmegaInCycles.*M>1))
        usFactorStr = sprintf('%0.3f,',deltaOmegaInCycles.*M);
        logMessage(['Warning: undersampling by a factor of (',usFactorStr,'). Image replication will occur! Reduce the number of pixels in the lateral dimension.']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = cztn(T,M,W,A,'centered');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The energy is spread over more sample points when sub-sampling. The next
    % lines makes sure that the input and output norm are the same when the output space size M goes to infinity.
    amplitudeCorrection = prod(sqrt(deltaOmegaInCycles(deltaOmegaInCycles~=0))); % But, do skip any singleton dimensions as the sample size would be undefined
%     subPixelOffsets = cellfun(@(outR) outR(1), outputRanges, 'UniformOutput',true);
%     amplitudeCorrection = amplitudeCorrection*exp(2i*pi*sum(subPixelOffsets));
    F = F*amplitudeCorrection;
    
    if nargout == 0
        plot(abs(F), 'Color', 'b'); hold on;
        plot(real(F), 'Color', 'g'); hold on;
        plot(imag(F), 'Color', 'r'); hold off;
        
        clear F;
    end
end
