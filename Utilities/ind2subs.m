%
% subs = ind2subs(sz, ind)
%
% Returns the same results as ind2sub but in a vector or matrix form.
%
% Inputs:
%     sz: a row vector indicating the size of the array
%     ind: one or more indices in the array
% 
% Output:
%     subs: a row vector of the same length as sz or a matrix consisting of
%     one such row per element of ind.
%
function subs = ind2subs(sz, ind)
    subs = zeros(numel(ind),numel(sz));
    if ~isempty(sz)
        nextIdx = ind(:);
        % loop invariant: [1 cumprod(sz(1:end-1))]*(subs'-1) + 1 == ind
        for dimIdx = 1:numel(sz)-1,
            currIdx = 1+mod(nextIdx-1,sz(dimIdx));
            nextIdx = 1+(nextIdx-currIdx)./sz(dimIdx);
            subs(:,dimIdx) = currIdx;
        end
        subs(:,end) = nextIdx;
    end
end