% shrunkVectors = shrink(l2Offsets, l1Weight, nbSphericalDimensions)
%
% Minimizes the sum of a weighted l1-norm and an l2-norm:
% l1Weight*sum(sqrt(sum(abs(x).^2,2))) + sum(abs(x(:)-l2Offsets(:)).^2) for x
%
% The soft thresholding algorithm [D.L. Donoho, De-noising by Soft-Thresholding, IEEE Trans. Inform. Theory, 41, 613–527, 1995]
% 
% l2Offsets = column vector indicating the L2 target, when multiple columns are
% specified, a vector variable is returned that minimizes the sum of the
% vector norms.
%
% nbSphericalDimensions = Optional: the number of dimensions to optimize jointly
% as sqrt(sum(abs(l2Offsets).^2,nbDimsOfOptimization+1)). Default: 1
%
% default: l1Weight=1
%
%
function shrunkVectors = shrink(l2Offsets, l1Weight, nbSphericalDimensions)
    if nargin<2,
        l1Weight=1;
    end
    if nargin<3 || isempty(nbSphericalDimensions),
        nbSphericalDimensions=1;
    end
    
    inputSize=size(l2Offsets);
    
    l2Offsets=reshape(l2Offsets,[],nbSphericalDimensions);
    
    if nbSphericalDimensions==1,
        R=abs(l2Offsets); % Just faster
    else
        R=sqrt(sum(abs(l2Offsets).^2,2));
    end
    % scaling always in the same quadrant
    scaling = max(0,1-l1Weight./(2*R+(R==0))); % When R==0, scaling doesn't matter
    
    shrunkVectors=l2Offsets.*repmat(scaling,[1 nbSphericalDimensions]);
    
    shrunkVectors=reshape(shrunkVectors,inputSize);
end