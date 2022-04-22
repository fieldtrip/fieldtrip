function [cfg] = ft_trialfun_balert(cfg)

% FT_TRIALFUN_BALERT extract trials from B-Alert data using an intermediate CSV file.
% FieldTrip cannot yet directly interpret the event markers from B-Alert data.
% Therefore, it is necessary to have B-Alert LAB. This is (paid) software from
% Advanced Brain Monitoring, in which you extract the eventmakers using the function:
% readevents(*.Events.edf, *.Signals.Raw.edf) to write a *.csv file.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset = string with the *.csv filename
%   cfg.trialfun = 'ft_trialfun_balert'
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

fid = fopen(cfg.dataset, 'rt');
frmt = '%d%d%s%d%d%s%d%d%d';
events = textscan(fid, frmt, 'Delimiter', ', ', 'CollectOutput', true, 'HeaderLines', 1);
len = size(events{1,1}, 1);
eventMarkers = str2num(cell2mat(events{1, 4}(2:len-1, :))); % event marker data
epochTime = events{1,1}(2:len-1, :); % epoched time
fclose(fid);

% Resample to 1024 Hz and save as sample numbers
cfg.event = struct(...
  'type', repmat('STATUS', size(epochTime, 1)), ...
  'sample', epochTime(:, 1)*1024 + round(epochTime(:, 2)/256*1024), ...
  'value', eventMarkers);

% Find stimulus
stimIndex = zeros(size(eventMarkers, 1), 1);
for i=1:length(cfg.eventvalue)
  stimIndex = stimIndex + (eventMarkers==cfg.eventvalue(i));
end
stimIndex = logical(stimIndex);
cfg.trl(:,1) = cfg.event.sample(stimIndex)-cfg.prestim*cfg.newfs;
cfg.trl(:,2) = cfg.event.sample(stimIndex)+cfg.poststim*cfg.newfs;
cfg.trl(:,3) = eventMarkers(stimIndex);


