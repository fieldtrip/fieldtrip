function [data] = avgref(data, sel)

% AVGREF computes the average reference in each column
%   [data] = avgref(data)
%
% or it computes the re-referenced data relative to the
% average over the selected channels
%   [data] = avgref(data, sel)

% Copyright (C) 1998-2002, Robert Oostenveld
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

% determine the dimension of the data
if length(size(data))==3
  % multiple epochs
  dim=3;
else
  % single epoch with multiple channels
  dim=2;
end

if nargin==1
  % default is to use all channels for average referencing
  if dim==3
    sel=1:size(data,2);
  else
    sel=1:size(data,1);
  end
end

if dim==3
  % the data contains multiple epochs
  for epoch=1:size(data,1)
    reference = mean(data(epoch,sel,:), 2);
    data(epoch,:,:) = data(epoch,:,:) - repmat(reference, [1 size(data,2) 1]);
  end
else
  % the data contains a single epoch
  reference = mean(data(sel,:), 1);
  data = data - repmat(reference, size(data,1), 1);
end

