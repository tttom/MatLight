 function data = poisson(xm, seed)
%function data = poisson(xm, seed)
%
% Generate Poisson random vector with mean xm.
% For small, use poisson1.m
% For large, use poisson2.m
% see num. rec. C, P. 222
%
% Copyright 1997-4-29, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && strcmp(xm, 'test'), poisson_test, return, end

if ~isvar('seed')
	seed = [];
end

data	= xm;
xm	= xm(:);

small = xm < 12;
data( small) = poisson1(xm( small), seed);
data(~small) = poisson2(xm(~small), seed);

%
% run a timing race against matlab's poissrnd
%
function poisson_test
n = 2^8;
t = reshape(linspace(1, 1000, n^2), n, n);
cpu tic
poisson(t);
cpu toc 'fessler poisson time'
tic 
poissrnd(t);
cpu toc 'matlab poissrnd time'

if 0 % look for bad values?
	poisson(10^6 * ones(n,n));
	poisson(1e-6 * ones(n,n));
	poisson(linspace(11,13,2^10+1));
	poisson(linspace(12700,91800,n^3));
end
