%
% XI = interpSinc(varargin)
%
% Interpolates a complex function on a regular grid using the CZT so the
% output ranges can be offset with respect to the input ranges.
% Input ranges can be left empty which defaults to [1:N], where N is the
% number of elements in that dimimension.
%
% Example: XI = interpSinc(xInputRange, yInputRange, X, xOutputRange, yOutputRange)
% The input size of X must be equal to the vector of the lengths of the
% preceeding input arguments. The size of the output XI is equal to the
% vector of the output ranges following the input argument X.
% (untested) X may have trailing dimensions that will be handled in
% parallel.
%
function XI = interpSinc(varargin)
    if nargin<1
%         inputFRange = 0.5+[-18:17].';
        inputFRange = exp(2) + [-18:5/pi:17].';
        X = max(0,inputFRange+10).*(inputFRange<15) + 10i*(abs(inputFRange)<10);

%         outputFRange = [-18:15].';
        dfRel = diff(inputFRange(1:2))*1.87;
        outputFRange = pi + [-18:dfRel:15].';
    
        varargin = {inputFRange, X, outputFRange};
    end

    % Parse inputs
    nbNumericArguments = numel(varargin);
    nbRanges = floor((nbNumericArguments-1)/2);
    arrayArgIdx = nbRanges+1;
    inputRanges = varargin(1:nbRanges);
    X = varargin{arrayArgIdx};
    outputRanges = varargin(arrayArgIdx+(1:nbRanges));
    
    % Fourier transform the input to the CZT
    if ndims(X) == nbRanges
        XIFT = ifftn(ifftshift(X));
    else
        for dim = 1:nbRanges
            XIFT = ifft(ifftshift(X,dim),[],dim);
        end
    end
    clear X;
    
    % Calculate the interpolated Fourier transform
    XI = interpFt(inputRanges{:}, XIFT, outputRanges{:});
   
    if nargin<1 && nargout<1
        figure;
        X = fftshift(fftn(XIFT));
        scatter(inputFRange, real(X), 'filled');
        hold on; scatter(inputFRange, imag(X), 'filled', 'FaceColor','r');
        hold on; plot(outputFRange, real(XI), 'Color','k', 'LineWidth', 2);
        hold on; plot(outputFRange, imag(XI), 'Color','r');
        legend({'orig real','orig imag','real','imag'}, 'Location','NorthWest');

        clear XI
    end
end
