function [cluster, numcluster] = findcluster(onoff, spatdimneighbstructmat, varargin)

% FINDCLUSTER returns all connected clusters for a three-dimensional six-connected
% neighborhood
%
% Use as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, minnbchan)
% or ar
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
% where
%   onoff                  =  is a 3D boolean matrix with size N1xN2xN3
%   spatdimneighbstructmat =  defines the neighbouring channels/combinations, see below
%   minnbchan              =  the minimum number of neighbouring channels/combinations
%   spatdimneighbselmat    =  is a special neighbourhood matrix that is used for selecting
%                             channels/combinations on the basis of the minnbchan criterium
%
% The neighbourhood structure for the first dimension is specified using
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column
% corresponds to a channel (combination) along the first dimension and along that
% row/column, elements with "1" define the neighbouring channel(s) (combinations).
% The first dimension of onoff should correspond to the channel(s) (combinations).
%
% See also SPM_BWLABEL, BWLABEL, BWLABELN

% Copyright (C) 2004-2020, Robert Oostenveld
% Copyright (C) 2021, Robert Oostenveld and Jan-Mathijs Schoffelen
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

% the calling code should ensure that SPM is on the path, preferably the latest version
ft_hastoolbox('spm', -1);

siz           = size(onoff);
spatdimlength = siz(1);
siz           = siz(2:end);
dims          = ndims(onoff);

% this is for reshaping back
celldims      = num2cell([siz 1]);

% this is for efficient matrix indexing
seldims = repmat({':'},[1 dims-1]);

% put the spatial dimension last: this substantially speeds up allocation
% and selection of data from the matrix, which outweighs the time spent to
% do the matrix permutation;
onoff = permute(onoff, [2:dims 1]);

if ~ismatrix(spatdimneighbstructmat) || ~all(size(spatdimneighbstructmat)==spatdimlength)
  ft_error('invalid dimension of spatdimneighbstructmat');
end

% this input argument handling is somewhat clunky, but it does the trick
if length(varargin)==1
  minnbchan = varargin{1};
else
  minnbchan = 0;
end
if length(varargin)==2
  spatdimneighbselmat = varargin{1};
  minnbchan           = varargin{2};
end

if minnbchan>0
  % For every (time,frequency)-element, it is calculated how many supra
  % threshold neighbours this spatial element has. If a supra threshold
  % spatial elementhas fewer than minnbchan supra threshold neighbours,
  % it is removed from onoff.
  
  if length(varargin)==1
    selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
  end
  if length(varargin)==2
    selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
  end
  
  nremoved = 1;
  while nremoved>0
    nsigneighb = reshape(reshape(single(onoff),[prod(siz) spatdimlength])*selectmat,[siz spatdimlength]);
    remove     = (onoff.*nsigneighb) < minnbchan;
    nremoved   = length(find(remove.*onoff));
    onoff(remove) = 0;
  end
end

% for each channel or channel-combination, find the connected time, frequency, or time-frequency clusters
labelmat = zeros(size(onoff));
numcluster = 0;

if ~(numel(siz)==1 && all(siz==1) && islogical(onoff))
  for j = 1:spatdimlength
    if numel(siz) <= 3 % if 2D or 3D data (without spatial dimension)
      % use SPM for 2D/3D data instead of the MATLAB image processing toolbox
      [clus, num] = spm_bwlabel(double(onoff(seldims{:},j)), 6);
    else
      [clus, num] = bwlabeln(double(onoff(seldims{:},j)), conndef(dims-1, 'min'));
    end
    clus(clus~=0) = clus(clus~=0) + numcluster;
    labelmat(seldims{:},j) = clus;
    numcluster = numcluster + num;
  end
else
  labelmat(onoff>0) = 1:sum(onoff(:));
  numcluster = sum(onoff(:));
end

% put the spatial dimension back upfront
labelmat = permute(labelmat, [dims 1:(dims-1)]);

% combine the non-spatial dimensions for simplicity
labelmat = reshape(labelmat, spatdimlength, []);

% combine clusters that are connected in neighbouring channels or channel
% combinations. Here we convert the input to uint32 as that is required by the mex
% file, and the values will be positive integers anyway.
if spatdimlength>1
  cluster = combineClusters(uint32(labelmat), logical(spatdimneighbstructmat), uint32(numcluster));
else
  cluster = labelmat;
end

% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, celldims{:});

% update the total number of clusters
numcluster = max(cluster(:));