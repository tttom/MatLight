% f = cztn(x, M, W, A, originCentered)
%
% Calculate the partial spectrum of x using the n-dimensional chirp z transform.
% x: input matrix, centered on central pixel. Implicit zero padding will
% occur at the right and will be corrected for.
% M, W, and A are vectors containing the scalars corresponding to each czt
% per dimension.
% originCentered: boolean, default false. Indicates if the input and output
% matrices have the origin in the central pixel (as defined by fftshift).
% Specify, true, 'centered' or 'originCentered' to set to true. Note that A
% remains the same.
%
% f: output matrix of size [M, size(x,length(M)+1) size(x,length(M)+2) ... size(x,ndims(x))]
%
function f = cztn(x, M, W, A, originCentered)
    if (nargin<2 || isempty(M))
        M = size(x);
    end
    if (nargin<3 || isempty(W))
        W = exp(-2i*pi./M);
    end
    if (nargin<4 || isempty(A))
        A = 1;
    end
    if (nargin<5 || isempty(originCentered))
        originCentered = false;
    end
    if (ischar(originCentered))
        switch(lower(originCentered))
            case {'centered','centred','origincentered','origincentred'}
                originCentered = true;
            otherwise
                originCentered = false;
        end     
    end

    nbDims = numel(M);
    
    % Work back to deltaRng
    maxFieldSpFreqInInputUnits = floor(size(x)/2);
    deltaRng = maxFieldSpFreqInInputUnits(1:nbDims).*imag(log(W))./(-2*pi);
    
    % Shift the output window by offsetting A
    if originCentered,
        %A = A.*W.^(floor(M/2)); % Not accurate enough!
        A = A.*exp(1i*angle(W).*floor(M./2)); % Only works for Fourier transforms
        % Adjust for z-transforms:
        zTransformDim = abs(abs(W)-1) > sqrt(eps(W));
        A(zTransformDim) = A(zTransformDim).*abs(W(zTransformDim)).^floor(M(zTransformDim)./2);
        
        outputPixelShift = angle(A)./angle(W);
    end
    
    f = x;
%     clear x;
    for dimIdx = 1:nbDims,
        if ~isnan(A(dimIdx)), % make sure it is not a dimension that has to be left alone
            if size(f,dimIdx)>1,
                f = cztw(f, M(dimIdx), W(dimIdx), A(dimIdx), dimIdx);
                if originCentered && deltaRng(dimIdx)~=0,
                    f = bsxfun(@times, f, shiftdim(exp(2i*pi*([1:size(f,dimIdx)]-1-outputPixelShift(dimIdx))*deltaRng(dimIdx)).',1-dimIdx)); %Correct the pre-shift induced phase error
                end
            else
                f = repmat(f,[ones(1,dimIdx-1) M(dimIdx) ones(1,ndims(f)-dimIdx)]);
            end
        end
    end
end
