% C = bsxVecFun(functor, A, B)
%
% Same as bsxfun but for functors that take vector inputs along the first
% dimension, such as @cross and @dot.
%
function C = bsxVecFun(functor, A, B)
    sizeA = size(A);
    sizeA(1) = 1;
    sizeB = size(B);
    sizeB(1) = 1;
    
    sizeA(end+1:numel(sizeB)) = 1;
    sizeB(end+1:numel(sizeA)) = 1;
    
    sizeC = max(sizeA,sizeB);
    
    A = repmat(A,sizeC./sizeA);
    B = repmat(B,sizeC./sizeB);
    
    C = functor(A,B);
end