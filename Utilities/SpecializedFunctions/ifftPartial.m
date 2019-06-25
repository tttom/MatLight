%
% T = ifftPartial(F, tRange, dim)
%
% Calculates the ifft for specific time points only
%
function T = ifftPartial(F, tRange, dim)
    if nargin<3 || isempty(dim),
        dim = 1;
    end
    T = conj(fftPartial(conj(F), tRange, dim))./size(F,dim);
end