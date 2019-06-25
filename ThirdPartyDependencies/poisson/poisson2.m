 function data = poisson2(xm, seed)
%function data = poisson2(xm, seed)
% Generate Poisson random column vector with mean xm.
% Uses rejection method - good for large values of xm.
% See "Numerical Recipes in C", P. 222.

if nargin < 2, help(mfilename), error(mfilename), end

if isvar('seed') & ~isempty(seed)
	rand('state', seed)
end

data = zeros(size(xm));
if any(xm < 0), error 'negative poisson means?', end
data(xm > 0) = poisson2_positive(xm(xm > 0));


%
% poisson2_positive()
%
function data = poisson2_positive(xm)
sx = sqrt(2.0 * xm);
lx = log(xm);
gx = xm .* lx - gammaln(1 + xm);

data = zeros(size(xm));
id = [1:length(xm)]'; % indicates which data left to do

%factor = 0.9; % from num rec
factor = 0.85; % seems to work better

while any(id)
	Tss = sx(id);
	Tll = lx(id);
	Tgg = gx(id);
	Txx = xm(id);

	yy = zeros(size(id));
	em = zeros(size(id));
	ib = true(size(id));

	while ib
		yy(ib) = tan(pi * rand(size(ib)));
		em(ib) = Tss(ib) .* yy(ib) + Txx(ib);
		ib = find(em < 0);
	end

	em = floor(em);
	tt = factor * (1+yy.*yy) .* exp(em .* Tll - gammaln(em+1) - Tgg);
	if any(tt > 1)
%		pr xm(tt > 1)
%		pr max(tt)
		error('factor is too large! please report xm value(s)')
	end

	ig = rand(size(id)) < tt;
	data(id(ig(:))) = em(ig);
	id = id(~ig(:));
end
