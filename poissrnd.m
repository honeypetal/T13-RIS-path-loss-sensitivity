function r = poissrnd(lambda, varargin)
% POISSRND  Poisson random numbers (stand-in for the Statistics Toolbox
% version, so SimRIS runs without that toolbox). Knuth's method.
%   r = poissrnd(lambda)           same size as lambda
%   r = poissrnd(lambda, m, n,...) lambda scalar, output m x n x ...
if nargin > 1
    sz = [varargin{:}];
    if isscalar(sz), sz = [sz sz]; end
    lambda = lambda * ones(sz);
end
r = zeros(size(lambda));
for k = 1:numel(lambda)
    L = exp(-lambda(k)); n = 0; p = rand;
    while p > L
        n = n + 1;
        p = p * rand;
    end
    r(k) = n;
end
end
