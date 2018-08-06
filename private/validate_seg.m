function [tissue1, tissue2, tissue3] = validate_seg(tissue1, tissue2, tissue3)

% VALIDATE_SEG ensures that the segmentation represents tissue types in a cumulative than exclusive 
% manner. 
%
% Use as
%   [tissue1, tissue2, tissue3] = validate_segmentation(tissue1, tissue2, tissue3)
% where the second two input (and output) arguments are optional. In case of more than one input 
% argument the tissue-types should follow eachother from inside towards outside (e.g. tissue1 = brain,
% tissue2 = skull, tissue = scalp). 
%
% The output will consist of one or more boolean segmentations without empty spaces inside. 
% In such way, more than one tissue-types will be represented in an overlapping manner. If
% the input is invalid and cannot be converted to overlapping segmentations, this function will give
% an error.
%
% This function makes use of functions from the MATLAB Signal Processing Toolbox.
%
% See also TRIANGULATE_SEG, PREPARE_MESH_SEGMENTATION

% Copyright (C) 2012, Robert Oostenveld
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

% the boundaries for surface extraction should be at
% t = sort(unique(brain+skull+2*scalp))
% t = (t(1:end-1) + t(2:end))/2

% FIXME ensure that there is no air in the brain


if false
  % this is an overlapping example
  tissue1 = [0 0 0 0 0 1 1 1 0 0 0 0];
  tissue2 = [0 0 0 1 1 1 1 1 1 0 0 0];
  tissue3 = [0 0 1 1 1 1 1 1 1 1 0 0];
  
  % this is an exclusive example
  tissue1 = [0 0 0 0 0 1 1 0 0 0 0 0];
  tissue2 = [0 0 0 1 1 0 0 1 1 0 0 0];
  tissue3 = [0 0 1 0 0 0 0 0 0 1 0 0];
  
  % this is an invalid/inconsistent example
  tissue1 = [0 0 0 1 0 0 0 0 0 0 0 0];
  tissue2 = [0 0 0 1 1 0 0 0 0 0 0 0];
  tissue3 = [0 0 0 0 0 0 1 0 0 0 0 0];
end


 if nargin<2
   % this default applies to an overlapping description
   tissue2 = tissue1;
 end

if nargin<3
  % this default applies to an overlapping description
  tissue3 = tissue2;
end



if ~isequal(size(tissue1), size(tissue2))
  ft_error('inconsistent size of segmentations')
  
elseif ~isequal(size(tissue1), size(tissue3))
  ft_error('inconsistent size of segmentations')
  
elseif ~isa(tissue1, 'logical') && ~all(tissue1(:)==0 | tissue1(:)==1)
  ft_error('the first tissue is not a binary segmentation');
  
elseif ~isa(tissue2, 'logical') && ~all(tissue2(:)==0 | tissue2(:)==1)
  ft_error('the second tissue is not a binary segmentation');
  
elseif ~isa(tissue3, 'logical') && ~all(tissue3(:)==0 | tissue3(:)==1)
  ft_error('the third tissue is not a binary segmentation');
  
end

% ensure that the first tissue is filled 
tissue1 = imfill(tissue1,'holes');  

if (~any(tissue1(:)&tissue2(:)))||(~any(tissue2(:)&tissue3(:)))||(~any(tissue3(:)&tissue1(:)))  
  % the segmentation is described as exclusive, i.e. there is no overlap
  % or it is a mixed representation
 

  % check for air inside the combination of tissue types
  air = ~imfill(tissue1|tissue2|tissue3, 'holes');
  
  if ~all(tissue1(:)|tissue2(:)|tissue3(:)|air(:))
    ft_error('there are voxels which do not belong to any tissue or air')
  end
  
  
  % convert them into an overlapping segmentation
  tissue1 = tissue1;
  tissue2 = tissue2 | tissue1;
  tissue3 = tissue3 | tissue2 | tissue1;


end

% the segmentation is described as overlapping
% this is suitable for detecting boundaries

if nargin > 1
  if any(tissue1(:)&~tissue2(:))
    ft_error('the first tissue outside of the second')
  end

  if any(tissue2(:)&~tissue3(:))
    ft_error('the second tissue is outside of the third')
  end

  if any(tissue1(:)&~tissue3(:))
    ft_error('there is first tissue is outside the third')
  end
  
  if ~any(tissue2(:)&~tissue1(:))
      ft_error('the first two tissues are not different')
  end
  
end

if nargin > 2
  if ~any(tissue3(:)&~tissue2(:))
      ft_error('the last two tissues are not different')
  end
end  

holes = imfill(tissue1, 'holes') & ~tissue1;
if any(holes(:))
  ft_error('there are holes in the first tissue');
end

if nargin > 1
 holes = imfill(tissue2, 'holes') & ~tissue2;
 if any(holes(:))
  ft_error('there are holes in the second tissue');
 end


 holes = imfill(tissue3, 'holes') & ~tissue3;
 if any(holes(:))
  ft_error('there are holes in the third tissue');
 end
end

