function [brain, skull, scalp] = validate_seg(brain, skull, scalp)

% VALIDATE_SEG ensures that the segmentation represents three tissue types
% that are overlapping rather than exclusive.
%
% Use as
%   [brain, skull, scalp] = validate_segmentation(brain, skull, scalp)
% where the skull and scalp input (and output) arguments are optional.
%
% The output will consist of three overlapping boolean segmentations. If
% the input is invalid and cannot be converted to overlapping
% segmentations, this function will give an error.
%
% See also TRIANGULATE_SEG

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% the boundaries for surface extraction should be at
% t = sort(unique(brain+skull+2*scalp))
% t = (t(1:end-1) + t(2:end))/2

% FIXME ensure that there is no air in the brain

if false
  % this is an overlapping example
  brain = [0 0 0 0 0 1 1 1 0 0 0 0];
  skull = [0 0 0 1 1 1 1 1 1 0 0 0];
  scalp = [0 0 1 1 1 1 1 1 1 1 0 0];
  
  % this is an exclusive example
  brain = [0 0 0 0 0 1 1 0 0 0 0 0];
  skull = [0 0 0 1 1 0 0 1 1 0 0 0];
  scalp = [0 0 1 0 0 0 0 0 0 1 0 0];
  
  % this is an invalid/inconsistent example
  brain = [0 0 0 1 0 0 0 0 0 0 0 0];
  skull = [0 0 0 1 1 0 0 0 0 0 0 0];
  scalp = [0 0 0 0 0 0 1 0 0 0 0 0];
end

if nargin<2
  % this default applies to an overlapping description
  skull = brain;
end

if nargin<3
  % this default applies to an overlapping description
  scalp = skull;
end


if ~isequal(size(brain), size(skull))
  error('inconsistent size of segmentations')
  
elseif ~isequal(size(brain), size(scalp))
  error('inconsistent size of segmentations')
  
elseif ~isa(brain, 'logical') && ~all(brain(:)==0 | brain(:)==1)
  error('the brain is not a binary segmentation');
  
elseif ~isa(skull, 'logical') && ~all(skull(:)==0 | skull(:)==1)
  error('the skull is not a binary segmentation');
  
elseif ~isa(scalp, 'logical') && ~all(scalp(:)==0 | scalp(:)==1)
  error('the scalp is not a binary segmentation');
  
elseif ~any(brain(:)&skull(:)) && ~any(skull(:)&scalp(:)) && ~any(scalp(:)&brain(:))
  % the segmentation is described as exclusive, i.e. there is no overlap
  
  % check for air inside the combination of tissue types
  air = ~imfill(brain|skull|scalp, 'holes');
  
  if ~all(brain(:)|skull(:)|scalp(:)|air(:))
    error('there are voxels which are not brain, skull, scalp or air')
  end
  
  % convert them into an overlapping segmentation
  brain = brain;
  skull = skull | brain;
  scalp = scalp | skull | brain;
  
end

% the segmentation is described as overlapping
% this is suitable for detecting boundaries

if any(brain(:)&~skull(:))
  error('there is brain outside the skull')
end

if any(skull(:)&~scalp(:))
  error('there is skull outside the scalp')
end

if any(brain(:)&~scalp(:))
  error('there is brain outside the scalp')
end

holes = imfill(brain, 'holes') & ~brain;
if any(holes(:))
  error('there are holes in the brain');
end

holes = imfill(skull, 'holes') & ~skull;
if any(holes(:))
  error('there are holes in the skull');
end

holes = imfill(scalp, 'holes') & ~scalp;
if any(holes(:))
  error('there are holes in the scalp');
end

