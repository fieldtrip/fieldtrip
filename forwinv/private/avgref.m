function [data] = avgref(data, sel);

% AVGREF computes the average reference in each column
%	[data] = avgref(data)
%
% or it computes the re-referenced data relative to the
% average over the selected channels
%	[data] = avgref(data, sel)

% Copyright (C) 1998-2002, Robert Oostenveld
%
% $Log: avgref.m,v $
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights
%

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
    reference = mean(squeeze(data(epoch,sel,:)), 1);
    data(epoch,:,:) = squeeze(data(epoch,:,:)) - repmat(reference, size(data,2), 1);
  end
else
  % the data contains a single epoch
  reference = mean(data(sel,:), 1);
  data = data - repmat(reference, size(data,1), 1);
end

