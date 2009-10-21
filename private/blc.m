function [data] = blc(data, interval, optional);

% BLC does a baseline correction using the prestimulus interval of the data
%
%   [data] = baseline(data, interval);
%   [data] = baseline(data, begin, end);
%
% If no begin and end are specified, the whole timeinterval is
% used for baseline correction.

% Copyright (C) 1998-2002, Robert Oostenveld 
%
% $Log: blc.m,v $
% Revision 1.3  2003/03/14 10:17:28  roberto
% fixed bug that was introduced by last change, changend from repmat to for-loop
%
% Revision 1.2  2003/03/13 16:44:45  roberto
% fixed bug with multiple epochs and single channel data
%

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

