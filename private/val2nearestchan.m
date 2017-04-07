function channame = val2nearestchan(data,val)

% VAL2NEARESTCHAN returns the label of the channel with the value nearest
% to the specified value.
%
% use as channame = val2nearestchan(data,val)
% val = [time y] with time in sec
% works only on raw data

% Copyright (C) 2009, Ingrid Nieuwenhuis
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

ft_defaults

[data] = ft_checkdata(data, 'datatype', 'raw');

% identify time point nearest to given time
timevec = [];
for iTr = 1:length(data.trial)
  timevec = [timevec data.time{iTr}];
end
  
tmp = nearest(timevec, val(1));
nearest_time = timevec(tmp);

for iTr = 1:length(data.trial)
  if ~isempty(data.time{iTr} == nearest_time)
    break
  end
end

trlNb = iTr;
sample = data.time{trlNb}==nearest_time;

% find nearest channel
chanNb = nearest(data.trial{trlNb}(:,sample),val(2));
channame = data.label{chanNb};


  
