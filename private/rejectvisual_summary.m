function [chansel, trlsel, cfg] = rejectvisual_summary(cfg, data);

% REJECTVISUAL_SUMMARY

% determine the initial selection of trials and channels
nchan = length(data.label);
ntrl  = length(data.trial);
cfg.channel = ft_channelselection(cfg.channel, data.label);
trlsel  = logical(ones(1,ntrl));
chansel = logical(zeros(1,nchan));
chansel(match_str(data.label, cfg.channel)) = 1;

% select the specified latency window from the data
% here it is done BEFORE filtering and metric computation
for i=1:ntrl
  begsample = nearest(data.time{i}, cfg.latency(1));
  endsample = nearest(data.time{i}, cfg.latency(2));
  data.time{i} = data.time{i}(begsample:endsample);
  data.trial{i} = data.trial{i}(:,begsample:endsample);
end

% compute the offset from the time axes
offset = zeros(ntrl,1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, data.fsample);
end

interactive = 1;
h = figure;
progress('init', cfg.feedback, 'computing metric');
level = zeros(sum(chansel), ntrl);
for i=1:ntrl
  progress(i/ntrl, 'computing metric %d of %d\n', i, ntrl);
  [dat, label, time, cfg.preproc] = preproc(data.trial{i}(chansel,:), data.label(chansel), data.fsample, cfg.preproc, offset(i));
  switch cfg.metric
    case 'var'
      level(:,i) = std(dat, [], 2).^2;
    case 'min'
      level(:,i) = min(dat, [], 2);
    case 'max'
      level(:,i) = max(dat, [], 2);
    case 'maxabs'
      level(:,i) = max(abs(dat), [], 2);
    case 'range'
      level(:,i) = max(dat, [], 2) - min(dat, [], 2);
    case 'kurtosis'
      level(:,i) = kurtosis(dat, [], 2);
    case '1/var'
      level(:,i) = 1./(std(dat, [], 2).^2);
    otherwise
      error('unsupported method');
  end
end
progress('close');
% reinsert the data for the selected channels
dum = nan*zeros(nchan, ntrl);
dum(chansel,:) = level;
level     = dum;
origlevel = dum;
clear dum

% start the interactive display of the data
while interactive
  % determine the maximum value
  level = origlevel;
  level(~chansel,:) = nan;
  level(:,~trlsel)  = nan;
  maxperchan = max(level,[],2);
  maxpertrl  = max(level,[],1);

  % make the three figures
  if gcf~=h, figure(h); end
  clf
  fprintf('selected %d trials from %d\n', sum(trlsel), ntrl);
  fprintf('selected %d channels from %d\n', sum(chansel), nchan);

  subplot(2,2,1);
  % imagesc(level(chansel, trlsel));
  imagesc(level);
  axis xy;
  % colorbar;
  title(cfg.method);
  ylabel('channel number');
  xlabel('trial number');

  subplot(2,2,2);
  plot(maxperchan,1:nchan, '.');
  abc = axis; axis([abc(1:2) 1 nchan]);
  ylabel('channel number');

  subplot(2,2,3);
  plot(1:ntrl, maxpertrl, '.');
  abc = axis; axis([1 ntrl abc(3:4)]);
  xlabel('trial number');

  [toggle, interactive] = smartinput('toggle the following trials [number, 0=interactive]: ', []);
  if interactive,
    if all(toggle==0)
      % use visual selection
      fprintf('make visual selection of trials to be removed...\n');
      [x, y] = select2d;
      toggle = find(1:ntrl>=x(1) & ...
        1:ntrl<=x(2) & ...
        maxpertrl(:)'>=y(1) & ...
        maxpertrl(:)'<=y(2));
    end
    trlsel(toggle) = ~trlsel(toggle);
    continue;
  end

  [toggle, interactive] = smartinput('toggle the following channels [number, 0=interactive]: ', []);
  if interactive,
    if all(toggle==0)
      % use visual selection
      fprintf('make visual selection of channels to be removed...\n');
      [x, y] = select2d;
      toggle = find(1:nchan>=y(1) & ...
        1:nchan<=y(2) & ...
        maxperchan(:)'>=x(1) & ...
        maxperchan(:)'<=x(2));
    end
    chansel(toggle) = ~chansel(toggle);
    continue;
  end

  [toggle, interactive] = smartinput('select the trials to be plotted seperately [number, 0=interactive]: ', []);
  if interactive,
    if all(toggle==0)
      % use visual selection
      fprintf('make visual selection of trials to be plotted seperately...\n');
      [x, y] = select2d;
      toggle = find(1:ntrl>=x(1) & ...
        1:ntrl<=x(2) & ...
        maxpertrl(:)'>=y(1) & ...
        maxpertrl(:)'<=y(2));
    end
    for i=1:length(toggle)
      figure
      % the data being displayed here is NOT filtered
      %plot(data.time{toggle(i)}, data.trial{toggle(i)}(chansel,:));
      plot(data.time{toggle(i)}, blc(data.trial{toggle(i)}(chansel,:)));
      title(sprintf('trial %d', toggle(i)));
    end
    continue;
  end

  [toggle, interactive] = smartinput('are you finished with your selection [Y, n]: ', 'y');
  if interactive,
    continue;
  end

end % while interactive

