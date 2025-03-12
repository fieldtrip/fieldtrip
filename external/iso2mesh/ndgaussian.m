function h = ndgaussian(r, sigma, ndim)
%
% h=ndgaussian(r, sigma, ndim)
%
% create an ND Gaussian or box filter kernel matrix
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    r: kernel half-width, the output is 2*r+1 in each dimension; if
%       missing, use 1
%    sigma: the standard deviation of the Gaussian; if not given, use 1; if
%       set to inf, output box filter
%    ndim: an integer for the output dimension; if not given, use 3
%
% output:
%    h: an ndim-dimensional matrix
%
% -- this function is part of the Iso2Mesh Toolbox (http://iso2mesh.sf.net)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if (nargin < 3)
    ndim = 3;
    if (nargin < 2)
        sigma = 1;
        if (nargin == 0)
            r = 1;
        end
    end
end
latt = cell(ndim, 1);
if (isinf(sigma)) % ND box filter
    if (ndim == 1)
        h = ones(1, (2 * r + 1));
    else
        h = ones((2 * r + 1) * ones(1, ndim));
    end
    return
end
[latt{:}] = meshgrid(-r:r);
latt = cellfun(@(x) x .* x, latt, 'uniformoutput', false);
h = exp(-(sum(cat(ndims(latt{1}) + 1, latt{:}), ndims(latt{1}) + 1)) ./ (2 * sigma * sigma));
h = h ./ (sum(h(:)));
