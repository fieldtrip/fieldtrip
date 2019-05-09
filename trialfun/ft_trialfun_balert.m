function [ cfg ] = ft_trialfun_balert( cfg )
%ft_definetrialBAlert: Extract trials B-Alert data from CSV file. 
%   FieldTrip cannot yet interpret the event markers from B-Alert data.
%   Therefore, it is necessary to have B-Alert LAB. This is (paid) software
%   from Advanced Brain Monitoring, in which you extract the evetmakers
%   using the function: readevents(*.Events.edf,*.Signals.Raw.edf) . This
%   function returns a *.csv file in which, for every epoch (second) 
    fid = fopen(cfg.dataset,'rt');
    frmt = '%d%d%s%d%d%s%d%d%d';
    events = textscan(fid,frmt,'Delimiter',',','CollectOutput',true, 'HeaderLines',1);
    len = size(events{1,1},1);
    eventMarkers = str2num(cell2mat(events{1,4}(2:len-1,:))); % event marker data
    epochTime = events{1,1}(2:len-1,:); %epoched time
    % Resample to 1024 Hz and save as sample numbers
    cfg.event = struct(...
        'type',repmat('STATUS',size(epochTime,1)),...
        'sample',epochTime(:,1)*1024 + round(epochTime(:,2)/256*1024),...
        'value',eventMarkers);

    % Find stimulus 
    stimIndex = zeros(size(eventMarkers,1),1);
    for i = 1:length(cfg.eventvalue)
        stimIndex = stimIndex + (eventMarkers==cfg.eventvalue(i));
    end
    stimIndex = logical(stimIndex);
    cfg.trl(:,1) = cfg.event.sample(stimIndex)-cfg.prestim*cfg.newfs;
    cfg.trl(:,2) = cfg.event.sample(stimIndex)+cfg.poststim*cfg.newfs;
    cfg.trl(:,3) = eventMarkers(stimIndex);
end

