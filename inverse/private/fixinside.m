function [source] = fixinside(source, target)

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
  target = 'logical';
end

if ~isfield(source, 'inside')
  if isfield(source, 'leadfield')
    % determine the positions outside the brain on basis of the empty cells
    source.inside = ~cellfun(@isempty, source.leadfield);
  elseif isfield(source, 'filter')
    % determine the positions outside the brain on basis of the empty cells
    source.inside = ~cellfun(@isempty, source.filter);
  elseif isfield(source, 'pos')
    % assume that all positions are inside the region of interest
    source.inside  = [1:size(source.pos,1)]';
    source.outside = [];
  elseif isfield(source, 'dim') && isfield(source, 'transform')
    % assume that all positions are inside the region of interest
    source.inside  = [1:prod(source.dim(1:3))]';
    source.outside = [];
  end
end

if ~isfield(source, 'inside')
  % nothing to do
  return;
end

% determine the format
if isa(source.inside, 'logical')
  current = 'logical';
elseif all(source.inside(:)==0 | source.inside(:)==1)
  source.inside = logical(source.inside);
  current = 'logical';
else
  current = 'indexed';
end


if strcmp(current, 'indexed') && strcmp(target, 'indexed')
  % nothing to do
elseif strcmp(current, 'logical') && strcmp(target, 'logical')
  % nothing to do
elseif strcmp(current, 'indexed') && strcmp(target, 'logical')
  % remove outside
  if isfield(source, 'outside')
    source = rmfield(source, 'outside');
  end
  % convert inside to a logical array
  if isfield(source, 'pos')
    tmp = false(size(source.pos,1), 1);
  elseif isfield(source, 'dim') && isfield(source, 'transform')
    tmp = false(source.dim(1:3));
  end
  tmp(source.inside) = true;
  source.inside = tmp;
elseif strcmp(current, 'logical') && strcmp(target, 'index')
  % convert to a vectors with indices
  tmp = source.inside;
  source.inside  = find( tmp(:));
  source.outside = find(~tmp(:));
else
  ft_error('incorrect specification of the insidestyle')
end

if strcmp(target, 'logical')
  if isfield(source, 'pos')
    % reshape it so that it matches the positions
    source.inside = reshape(source.inside, size(source.pos,1), 1);
  elseif isfield(source, 'dim') && isfield(source, 'transform')
    % reshape it so that it matches the volume dimensions
    source.inside = reshape(source.inside, source.dim(1:3));
  end
end
