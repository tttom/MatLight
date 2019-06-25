%
% FI = interpFt(varargin)
%
% Calculates the interpolated Fourier Transform of a complex function on a
% regular grid using the CZT so that the output ranges can be offset with
% respect to the input ranges.
% Input ranges can be left empty which defaults to [1:N], where N is the
% number of elements in that dimimension.
% T should be the Inverse Fourier transform of the ifftshited to-be-interpolated.
% input, i.e. T = ifftn(ifftshift(F))
%
% Example: FI = interpFt(xInputRange, yInputRange, T, xOutputRange, yOutputRange)
% The input size of T must be equal to the vector of the lengths of the
% preceeding input arguments. The size of the output FI is equal to the
% vector of the output ranges following the input argument T.
% (untested) T may have trailing dimensions that will be handled in
% parallel.
%
function FI = interpFt(varargin)
    if nargin<1,
        inputFRange = exp(2)+[-18:5/pi:17].';
        F = max(0,inputFRange+10).*(inputFRange<15) + 10i*(abs(inputFRange)<10);
        T = ifftn(ifftshift(F));

        dfRel = diff(inputFRange(1:2))*.01;
        outputFRange = pi+[-15:dfRel:14].';
    
        varargin = {inputFRange, T, outputFRange};
    end

    nbNumericArguments = numel(varargin);
    nbRanges = floor((nbNumericArguments-1)/2);
    arrayArgIdx = nbRanges+1;
    inputRanges = varargin(1:nbRanges);
    T = varargin{arrayArgIdx};
    outputRanges = varargin(arrayArgIdx+(1:nbRanges));
    
    % Set empty ranges to default value
    defaultInputRanges = arrayfun(@(len) (1:len), size(T), 'UniformOutput',false);
    function rng = fixRange(rng, defaultRng)
        if isempty(rng),
            rng = defaultRng;
        end
    end
    inputRanges = cellfun(@fixRange, inputRanges, defaultInputRanges(1:numel(inputRanges)), 'UniformOutput',false);
    outputRanges = cellfun(@fixRange, outputRanges, inputRanges, 'UniformOutput',false);

    % Center DC component for using the CZT for interpolation
    if ndims(T) == nbRanges,
        T = fftshift(T);
    else
        for dim = 1:nbRanges,
            T = fftshift(T, dim);
        end
    end
    
    % Calculate relative outputRanges for the CZT input arguments
    function d = getD(vec)
        if numel(vec)>=2,
            d = diff(vec(1:2));
        else
            d = 1;
        end
    end
    relOutputRanges = cellfun(@(inR, outR) (outR - inR(1+floor(end/2)))./getD(inR),...
        inputRanges, outputRanges, 'UniformOutput',false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FI = cztFromRanges(T, relOutputRanges{:});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % correct scaling
    inputRangeLengths = cellfun(@(rIn) numel(rIn), inputRanges, 'UniformOutput',true);
    magnification = cellfun(@(rIn, rOut) getD(rOut)/getD(rIn),...
        inputRanges, outputRanges, 'UniformOutput',true);
    
    FI = FI * sqrt(prod(inputRangeLengths)/prod(magnification));
    
    if nargin<1 && nargout<1,
        figure;
        F = fftshift(fftn(ifftshift(T)));
        scatter(inputFRange, real(F), 'filled', 'FaceColor','k');
        hold on; scatter(inputFRange, imag(F), 'filled', 'FaceColor','r');
        hold on; plot(outputFRange, real(FI), 'Color','k', 'LineWidth', 2);
        hold on; plot(outputFRange, imag(FI), 'Color','r');
        legend({'orig real','orig imag','real','imag'}, 'Location','NorthWest');

        clear FI
    end
end
