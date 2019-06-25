%
% F = fftPartial(T, fRange, dim)
%
% Applies the 1D Fourier transform at a few specific frequencies.
%
function F = fftPartial(T, fRange, dim)
    if nargin<3 || isempty(dim),
        dim = 1;
    end

    if dim==1,
        inputSize = size(T);
        tRange = ifftshift([1:inputSize(1)]-1-floor(inputSize(1)/2))./inputSize(1);
        FT = exp(-2i*pi*fRange(:)*tRange);
        F = FT*T(:,:);
        F = reshape(F, [numel(fRange), inputSize(2:end)]);
    else
        % flip dimensions, delegate, and reorder
        nbDims = max(dim, ndims(T));
        P = [dim, 1:dim-1, dim+1:nbDims];
        T = permute(T, P);
        F = fftPartial(T, fRange);
        F = ipermute(F, P);
    end
end