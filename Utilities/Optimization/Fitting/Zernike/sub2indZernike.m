% js=sub2indZernike(n,m)
%
% Converts a Zernike coordinate (n,m), js, to an indexes 
% When multiple values are specified, js will have the same shape as
% the inputs n and m.
%
% The ordering is described here:
% Noll, R. J. (1976). "Zernike polynomials and atmospheric turbulence" (PDF). J. Opt. Soc. Am. 66 (3): 207. Bibcode:1976JOSA...66..207N. doi:10.1364/JOSA.66.000207.
%
% See also: ind2subZernike
%
function js=sub2indZernike(n,m)
    if nargin<2 || isemtpy(m),
        m=0;
    end
    js = n.*(n+1)./2; % number up to n-1
    js = js + abs(m) + (m==0); % correct number or one too low
    js = js + double(m~=0&xor(m<0,mod(js,2))); % make js odd if m negative
end
