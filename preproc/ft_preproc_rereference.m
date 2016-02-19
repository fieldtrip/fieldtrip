function [dat, ref] = ft_preproc_rereference(dat, refchan, method)

% FT_PREPROC_REREFERENCE computes the average reference over all EEG channels
% or rereferences the data to the selected channels
%
% Use as
%   [dat] = ft_preproc_rereference(dat, refchan)
% where
%   dat        data matrix (Nchans X Ntime)
%   refchan    vector with indices of the new reference channels
%   method     string, can be 'avg' or 'median'
%
% If the new reference channel is not specified, the data will be
% rereferenced to the average of all channels.
%
% See also PREPROC

% Copyright (C) 1998-2008, Robert Oostenveld
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

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the new reference channels
if nargin<2 || isempty(refchan) || (ischar(refchan) && strcmp(refchan, 'all'))
  refchan = 1:Nchans;
end
if nargin<3 || isempty(method)
  method = 'avg';
end

% compute the average value over the reference channels
if strcmp(method, 'avg')
    ref = mean(dat(refchan,:), 1);
elseif strcmp(method, 'median')
    ref = median(dat(refchan,:), 1);
end

% apply the new reference to the data
for chan=1:Nchans
  dat(chan,:) = dat(chan,:) - ref;
end
