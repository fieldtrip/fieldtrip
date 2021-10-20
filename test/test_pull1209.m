function test_pull1209

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY xdf2fieldtrip

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/pull1209.xdf');

% read the data and events
[data, event] = xdf2fieldtrip(filename);

% pull the EEG signals towards the baseline
data.trial{1} = ft_preproc_baselinecorrect(data.trial{1});

% start the time axis at t=0
data.time{1} = data.time{1} - data.time{1}(1);

%%

figure
plot(data.time{1}, data.trial{1});
ax = axis;

%%

for i=1:numel(event)
  % compute the sample number corresponding to the LSL timestamp
  event(i).sample = round((event(i).timestamp - data.hdr.FirstTimeStamp)/data.hdr.TimeStampPerSample) + 1;
  
  try
    x(1) = data.time{1}(event(i).sample);
    x(2) = data.time{1}(event(i).sample);
    y(1) = ax(3);
    y(2) = ax(4);
    
    h = line(x, y);
    set(h, 'Color', 'r')
  catch
    % events can fall outside the range of the continuous data
  end
  
end
