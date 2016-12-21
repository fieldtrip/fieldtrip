function [source] = fixinside(source, opt)

% FIXINSIDE ensures that the region of interest (which is indicated by the
% field "inside") is consistently defined for source structures and volume
% structures. Furthermore, it solves backward compatibility problems.
%
% Use as
%   [source] = fixinside(source, 'logical');
% or
%   [source] = fixinside(source, 'index');

% Copyright (C) 2006, Robert Oostenveld
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


if nargin<2
  opt = 'logical';
end

if ~isfield(source, 'inside')
  if isfield(source, 'pos')
    % assume that all positions are inside the region of interest
    source.inside  = [1:size(source.pos,1)]';
    source.outside = [];
  elseif isfield(source, 'dim')
    source.inside  = [1:prod(source.dim)]';
    source.outside = [];
  end
end

if ~isfield(source, 'inside')
  % nothing to do
  return;
end

% determine the format
if isa(source.inside, 'logical')
  logicalfmt = 1;
elseif all(source.inside(:)==0 | source.inside(:)==1)
  source.inside = logical(source.inside);
  logicalfmt = 1;
else
  logicalfmt = 0;
end

if ~logicalfmt && strcmp(opt, 'logical')
  % convert to a logical array
  if ~isfield(source, 'outside')
    source.outside = [];
  end
  if isfield(source, 'pos')
    tmp  = false(size(source.pos,1),1);
  elseif isfield(source, 'dim')
    tmp  = false(prod(source.dim),1);
  end
  tmp(source.inside) = true;
  if isfield(source, 'outside')
    tmp(source.outside) = false;
    source = rmfield(source, 'outside');
  end
  source.inside = tmp(:);
elseif logicalfmt && strcmp(opt, 'index')
  % convert to a vectors with indices
  tmp = source.inside;
  source.inside  = find( tmp(:));
  source.outside = find(~tmp(:));
else
  % nothing to do
end
