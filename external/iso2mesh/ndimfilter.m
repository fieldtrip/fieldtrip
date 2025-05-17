function img = ndimfilter(im, kernel, varargin)
%
% img=ndimfilter(im,kernel,r,sigma)
%
% filter an ND array using a specified filter using convolution
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    im: input ND array
%    kernel: can be an ND array, or a string. if string, the below filters
%        are supported:
%        'box': box filter (need r)
%        'gaussian': Gaussian filter (need r,sigma input)
%    r: kernel half-width, the output is 2*r+1 in each dimension; if
%       missing, use 1
%    sigma: the standard deviation of the Gaussian; if not given, use 1; if
%       set to inf, output box filter
%
% output:
%    img: the filtered ND array
%
% -- this function is part of the Iso2Mesh Toolbox (http://iso2mesh.sf.net)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if (nargin < 2)
    kernel = 'box';
end
if (ischar(kernel))
    switch (kernel)
        case 'box'
            kernel = ndgaussian(varargin{1}, inf, ndims(im));
        case 'gaussian'
            kernel = ndgaussian(varargin{1}, varargin{2}, ndims(im));
        otherwise
            error('filter type %s is not supported', type);
    end
end

img = convn(im, kernel, 'same');
