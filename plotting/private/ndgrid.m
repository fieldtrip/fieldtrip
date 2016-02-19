function [varargout] = ndgrid(varargin)

% NDGRID Generation of arrays for N-D functions and interpolation.
% [X1,X2,X3,...] = NDGRID(x1,x2,x3,...) transforms the domain
% specified by vectors x1,x2,x3, etc. into arrays X1,X2,X3, etc. that
% can be used for the evaluation of functions of N variables and N-D
% interpolation.  The i-th dimension of the output array Xi are copies
% of elements of the vector xi.
%
% [X1,X2,...] = NDGRID(x) is the same as [X1,X2,...] = NDGRID(x,x,...).
%
% For example, to evaluate the function  x2*exp(-x1^2-x2^2-x^3) over the
% range  -2 < x1 < 2,  -2 < x2 < 2, -2 < x3 < 2,
%
%     [x1,x2,x3] = ndgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%     z = x2 .* exp(-x1.^2 - x2.^2 - x3.^2);
%     slice(x2,x1,x3,z,[-1.2 .8 2],2,[-2 -.2])
%
% NDGRID is like MESHGRID except that the order of the first two input
% arguments are switched (i.e., [X1,X2,X3] = NDGRID(x1,x2,x3) produces
% the same result as [X2,X1,X3] = MESHGRID(x2,x1,x3)).  Because of
% this, NDGRID is better suited to N-D problems that aren't spatially
% based, while MESHGRID is better suited to problems in cartesian
% space (2-D or 3-D).
%
% This is a drop-in replacement for the MATLAB version in elmat, which is
% relatively slow for big grids. Note that this function only works up
% to 5 dimensions
%
% See also MESHGRID, INTERPN.

% Copyright(C) 2010, Jan-Mathijs Schoffelen, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin==0
  error('MATLAB:ndgrid:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin==1, varargin = repmat(varargin,[1 max(nargout,2)]); end

ndims = numel(varargin);
switch ndims
case 2
  ones1 = ones(1,numel(varargin{1}));
  ones2 = ones(1,numel(varargin{2}));
  
  x   = varargin{1}(:);
  y   = varargin{2}(:)';
  
  varargout{1} = x(:, ones2);
  varargout{2} = y(ones1, :);
case 3
  ones1 = ones(1,numel(varargin{1}));
  ones2 = ones(1,numel(varargin{2}));
  ones3 = ones(1,numel(varargin{3}));
  
  x   = varargin{1}(:);
  y   = varargin{2}(:)';
  z   = zeros(1,1,numel(varargin{3}));
  z(:) = varargin{3};
  
  varargout{1} = x(:, ones2, ones3);
  varargout{2} = y(ones1, :, ones3);
  varargout{3} = z(ones1, ones2, :);
case 4
  ones1 = ones(1,numel(varargin{1}));
  ones2 = ones(1,numel(varargin{2}));
  ones3 = ones(1,numel(varargin{3}));
  ones4 = ones(1,numel(varargin{4}));
  
  x   = varargin{1}(:);
  y   = varargin{2}(:)';
  z   = zeros(1,1,numel(varargin{3}));
  z(:) = varargin{3};
  xx   = zeros(1,1,1,numel(varargin{4}));
  xx(:) = varargin{4};
  
  varargout{1} = x(:, ones2, ones3, ones4);
  varargout{2} = y(ones1, :, ones3, ones4);
  varargout{3} = z(ones1, ones2, :, ones4);
  varargout{4} = xx(ones1, ones2, ones3, :);
case 5
  ones1 = ones(1,numel(varargin{1}));
  ones2 = ones(1,numel(varargin{2}));
  ones3 = ones(1,numel(varargin{3}));
  ones4 = ones(1,numel(varargin{4}));
  ones5 = ones(1,numel(varargin{5}));
  
  x   = varargin{1}(:);
  y   = varargin{2}(:)';
  z   = zeros(1,1,numel(varargin{3}));
  z(:) = varargin{3};
  xx   = zeros(1,1,1,numel(varargin{4}));
  xx(:) = varargin{4};
  yy   = zeros(1,1,1,1,numel(varargin{5}));
  yy(:) = varargin{5};
  
  varargout{1} = x(:, ones2, ones3, ones4, ones5);
  varargout{2} = y(ones1, :, ones3, ones4, ones5);
  varargout{3} = z(ones1, ones2, :, ones4, ones5);
  varargout{4} = xx(ones1, ones2, ones3, :,ones5);
  varargout{5} = yy(ones1, ones2, ones3, :,ones5);
otherwise
  error('this version of ndgrid supports inputs up to 5 dimensions');
  %call the ndgrid from elmat
  %FIXME this has to be done
end
