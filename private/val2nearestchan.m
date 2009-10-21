function channame = val2nearestchan(data,val)

% VAL2NEARESTCHAN returns the label of the channel with the value nearest
% to the specified value.
%
% use as channame = val2nearestchan(data,val)
% val = [time y] with time in sec
% works only on raw data

% Copyright (C) 2009, Ingrid Nieuwenhuis
%
% $Log: val2nearestchan.m,v $
% Revision 1.1  2009/10/09 12:05:06  ingnie
% first implementation, used by databrowser
%

fieldtripdefs

[data] = checkdata(data, 'datatype', 'raw');

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


  