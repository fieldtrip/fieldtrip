function [data] = blc(data, interval, optional)

% BLC does a baseline correction using the prestimulus interval of the data
%
%   [data] = baseline(data, interval);
%   [data] = baseline(data, begin, end);
%
% If no begin and end are specified, the whole timeinterval is
% used for baseline correction.

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
  nepoch = size(data,1);
  nchan  = size(data,2);
  ntime  = size(data,3);
else
  % single epoch with multiple channels
  dim=2;
  nchan  = size(data,1);
  ntime  = size(data,2);
end

% determine the interval to use for baseline correction
if nargin==1
  % default is to use the whole timeinterval for baseline correction
  interval = 1:ntime;
elseif nargin==2
  % use the interval specified by the user
elseif nargin==3
  % create the interval from the specified begin and end point
  interval = interval:optional;
end

if dim==3
  % the data contains multiple epochs
  for epoch=1:nepoch
    baseline = mean(data(epoch,:,interval), 3);
    for chan=1:nchan
      data(epoch,chan,:) = data(epoch,chan,:) - baseline(chan);
    end
  end
else
  % the data contains a single epoch
  baseline = mean(data(:,interval), 2);
  for chan=1:nchan
    data(chan,:) = data(chan,:) - baseline(chan);
  end
end

