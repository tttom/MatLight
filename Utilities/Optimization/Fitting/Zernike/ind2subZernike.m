% [n,m] = ind2subZernike(js)
%
% Converts a Zernike index, js>0, to a pair of subscript m,n, 0 <= m <= n
% When multiple values are specified, m and n will have the same shape as
% the input js.
%
% The ordering is described here:
% Noll, R. J. (1976). "Zernike polynomials and atmospheric turbulence" (PDF). J. Opt. Soc. Am. 66 (3): 207. Bibcode:1976JOSA...66..207N. doi:10.1364/JOSA.66.000207.
%
% See also: sub2indZernike
%
function [n,m] = ind2subZernike(js)
%     n=floor((1+sqrt(8*js-7))./2)-1;
    n=ceil((sqrt(1+8*js)-1)/2)-1;
    if nargout>1,
        mSeq=js-n.*(n+1)./2-1; % the zero-based sequence number for the real m = 0, 2, -2, 4, -4, 6, -6,... or 1, -1, 3, -3, 6, -6, ... (or the inverse depending on mod(j,2) )
        m=2*floor((mSeq+(1-mod(n,2)))./2)+mod(n,2); % absolute value of real m
        m=m.*(1-2*mod(js,2)); % If j odd, make m negative.
    end
end
