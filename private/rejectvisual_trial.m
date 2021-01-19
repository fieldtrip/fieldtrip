function [chansel, trlsel, cfg] = rejectvisual_trial(cfg, data)

% SUBFUNCTION for ft_rejectvisual

% determine the initial selection of trials
ntrl = length(data.trial);
if isequal(cfg.trials, 'all') % support specification like 'all'
  cfg.trials = 1:ntrl;
elseif isempty(cfg.trials)
  cfg.trials = 1:ntrl;
elseif isnumeric(cfg.trials)
  % use the selection as it is
end
trlsel = false(1, ntrl);
trlsel(cfg.trials) = true;

% determine the initial selection of channels
nchan = length(data.label);
cfg.channel = ft_channelselection(cfg.channel, data.label); % support specification like 'all'
chansel = false(1, nchan);
chansel(match_str(data.label, cfg.channel)) = true;

% compute the sampling frequency from the first two timepoints
fsample = 1/mean(diff(data.time{1}));

% compute the offset from the time axes
offset = zeros(ntrl,1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, fsample);
end

if (isfield(cfg, 'preproc') && ~isempty(cfg.preproc))
  ft_progress('init', cfg.feedback, 'filtering data');
  for i=1:ntrl
    ft_progress(i/ntrl, 'filtering data in trial %d of %d\n', i, ntrl);
    [data.trial{i}, label, time, cfg.preproc] = preproc(data.trial{i}, data.label, data.time{i}, cfg.preproc);
  end
  ft_progress('close');
end

% select the specified latency window from the data
% this is done AFTER the filtering to prevent edge artifacts
if ischar(cfg.latency)
  cfg.latency(1) = min(cellfun(@min, data.time));
  cfg.latency(2) = max(cellfun(@max, data.time));
end

% the info structure will be attached to the figure
% and passed around between the callback functions
info.trlop   = 1;
info.chanlop = 1;
info.quit    = 0;
info.ntrl    = ntrl;
info.nchan   = nchan;
info.data    = data;
info.cfg     = cfg;
info.offset  = offset;
info.chansel = chansel;
info.trlsel  = trlsel;

tmpcfg = [];
tmpcfg.layout = 'ordered';
tmpcfg.skipscale = 'yes';
tmpcfg.skipcomnt = 'yes';
for i=1:nchan
  tmpcfg.channel{i} = sprintf('channel %s', data.label{i});
end
info.layout = ft_prepare_layout(tmpcfg);

minx = min(info.layout.pos(:,1) - info.layout.width/2);
maxx = max(info.layout.pos(:,1) + info.layout.width/2);
miny = min(info.layout.pos(:,2) - info.layout.height/2);
maxy = max(info.layout.pos(:,2) + info.layout.height/2);

h = create_figure([minx maxx miny maxy], info);

% Maintain focus while GUI is open
while ishandle(h)
  drawnow
  info = guidata(h);
  if info.quit == 0
    uiwait;
  else
    chansel = info.chansel;
    trlsel = info.trlsel;
    delete(h);
    break
  end
end % while ishandle(h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function markexclude_next(h, event)
markexclude(h, event);
step_trial(h, event, 1);
end

function markinclude_next(h, event)
markinclude(h, event);
step_trial(h, event, 1);
end

function step_trial(h, event, step)
info = guidata(h);
if info.trlop == 1 && step < 0
  fprintf('At first trial\n');
elseif info.trlop == info.ntrl && step > 0
  fprintf('At last trial\n');
else
  trls = (1:info.ntrl);
  info.trlop = trls(nearest(trls, info.trlop + step));  
  % Update plots
  new_plots(info);
end
guidata(h, info);
uiresume;
end

function markinclude(h, event)
info = guidata(h);
% Include trial if it's not already included
if ~info.trlsel(info.trlop)
  info.trlsel(info.trlop) = 1;
  fprintf(description_trial(info));
  title(description_trial(info));
  % Color in black all channels that are included (excluded remain blank)
  set(info.h_chan(info.chansel), 'Color', 'k');
  guidata(h,info);
end
uiresume;
end

function markexclude(h, event)
info = guidata(h);
% Exclude channel if it's not already excluded
if info.trlsel(info.trlop)
  info.trlsel(info.trlop) = 0;
  fprintf(description_trial(info));
  title(description_trial(info));
  % Blank all channels
  set(info.h_chan(info.chansel), 'Color', 'none');
  guidata(h,info);
end
uiresume;
end

function key(h, event)

switch lower(event.Key)
  case 'rightarrow'
    step_trial(h, event, 1);
  case 'leftarrow'
    step_trial(h, event, -1);
  case 'g'
    markinclude(h, event);
  case 'b'
    markexclude(h, event);
  case 'q'
    stop(h);
  otherwise
    fprintf('unknown key pressed\n');
end
end

function button(h, event)
% Find selected trial
pos = get(gca, 'CurrentPoint');
x = pos(1,1);
y = pos(1,2);
info = guidata(h);
dx = info.layout.pos(:,1) - x;
dy = info.layout.pos(:,2) - y;
dd = sqrt(dx.^2 + dy.^2);
[d, i] = min(dd);
threshold = median(sqrt(info.layout.width.^2 + info.layout.height.^2)/2);
if d<threshold && i<=info.nchan
  % Toggle
  info.chansel(i) = ~info.chansel(i); % toggle
  info.chanlop = i;
  guidata(h,info);
    
  % Print to command window
  fprintf(description_channel(info));
  
  % Update plots
  update_plots(info);
  
  % Update trial color
  if strcmp(info.h_chan(i).Color, 'none') % channel is toggled to include
    info.h_chan(i).Color = 'k';
  else % channel is toggled to exclude
    info.h_chan(i).Color = 'none';
  end
else
  fprintf('button clicked\n');
end
uiresume;
end

function stop(h, event)
info = guidata(h);
info.quit = 1;
guidata(h,info);
uiresume;
end

function str = description_channel(info)
if info.chansel(info.chanlop)
  str = sprintf('channel %s marked to INCLUDE\n', info.data.label{info.chanlop});
else
  str = sprintf('channel %s marked to EXCLUDE\n', info.data.label{info.chanlop});
end
end

function str = description_trial(info)
if info.trlsel(info.trlop)
  str = sprintf('trial %d marked to INCLUDE\n', info.trlop);
else
  str = sprintf('trial %d marked to EXCLUDE\n', info.trlop);
end
end

function h = create_figure(limits, info)
% Creates the GUI and plots the data for the first time.
% All additional changes are done by manipulating the properties of the
% specific graphical objects
%
% Input:
% limits - limits of the area where the trials are plotted
% info   - the struct that is attched to figure using guidata

% Create figure
h = figure();
axis('off')
axis(limits);

% Set buttons and callback functions
info.ui.quit        = uicontrol(h,'units','pixels','position',[  5 5 40 18],'String','quit','Callback',@stop);
info.ui.prev        = uicontrol(h,'units','pixels','position',[ 50 5 25 18],'String','<','Callback',{@step_trial, -1});
info.ui.next        = uicontrol(h,'units','pixels','position',[ 75 5 25 18],'String','>','Callback',{@step_trial, 1});
info.ui.prev10      = uicontrol(h,'units','pixels','position',[105 5 25 18],'String','<<','Callback',{@step_trial, -10});
info.ui.next10      = uicontrol(h,'units','pixels','position',[130 5 25 18],'String','>>','Callback',{@step_trial, 10});
info.ui.exclude     = uicontrol(h,'units','pixels','position',[160 5 70 18],'String','exclude','Callback',@markexclude);
info.ui.include     = uicontrol(h,'units','pixels','position',[230 5 70 18],'String','include','Callback',@markinclude);
info.ui.excludenext = uicontrol(h,'units','pixels','position',[310 5 70 18],'String','exclude >','Callback',@markexclude_next);
info.ui.includenext = uicontrol(h,'units','pixels','position',[380 5 70 18],'String','include >','Callback',@markinclude_next);
set(h, 'WindowButtonUpFcn', @button);
set(h, 'KeyPressFcn', @key);

% Determine the y-axis limits for this channel for all included trials
[ymin, ymax] = get_ylim(info);

% Plot individual channels
hold('on')
color = 'k';
tim = info.data.time{info.trlop};
for chanindx=1:info.nchan  
  dat = info.data.trial{info.trlop}(chanindx,:);
  info.h_chan(chanindx) = ft_plot_vector(tim, dat, 'hpos', info.layout.pos(chanindx,1), 'vpos', info.layout.pos(chanindx,2), 'width', info.layout.width(chanindx), 'height', info.layout.height(chanindx), 'vlim', [ymin ymax], 'box', istrue(info.cfg.box), 'color', color, 'label', info.layout.label{chanindx}, 'labelpos', 'lowerleft');
end

% Only blank excluded channels if the trial is included
if info.trlsel(info.trlop)
  set(info.h_chan(~info.chansel), 'Color', 'none');
else % If channel is excluded blank out everything
  set(info.h_chan, 'Color', 'none');
end
hold('off')

% Enable or disable buttons as appropriate
check_buttons_status(info);

% Add title
title(description_trial(info),'interpreter','none');

% Update the GUI
guidata(h, info);
end

function [ymin, ymax] = get_ylim(info)
% Determines the y-axis limits for one trial for all included channels

if ~isnumeric(info.cfg.ylim)
  ymin = nan(1,info.nchan);
  ymax = nan(1,info.nchan);
  for chanindx = find(info.chansel)
    dat = info.data.trial{info.trlop}(chanindx,:);
    ymin(chanindx) = min(dat);
    ymax(chanindx) = max(dat);
  end
  ymin = min(ymin);
  ymax = max(ymax);
  if strcmp(info.cfg.ylim, 'maxabs') % handle maxabs, make y-axis center on 0
    ymax = max(abs(ymax), abs(ymin));
    ymin = -ymax;
  elseif strcmp(info.cfg.ylim, 'zeromax')
    ymin = 0;
  elseif strcmp(info.cfg.ylim, 'minzero')
    ymax = 0;
  end
else
  ymin = info.cfg.ylim(1);
  ymax = info.cfg.ylim(2);
end
end

function update_plots(info)
% Updates the trial plots after inclusion/exclusion of trials
%
% Input:
% info     - info struct with all the data attched to the GUI

[ymin, ymax] = get_ylim(info);

% Get positions and heights
v_centers = info.layout.pos(:, 2);
heights   = info.layout.height;

for chanindx = find(info.chansel)
  % Shift the vertical axis to zero
  vdat = info.data.trial{info.trlop}(chanindx, :); 
  % Scale to length 1 of the new data
  vdat = vdat ./ (ymax - ymin);
  % Scale by availabe plot area
  vdat = vdat .* heights(chanindx);
  % Shift to the vertical position
  vdat = vdat + v_centers(chanindx);
  % Update
  info.h_chan(chanindx).YData = vdat;
end
end

function new_plots(info)
% Updates the trial plots after channel change
%
% Input:
% info - info struct with all the data attched to the GUI

[ymin, ymax] = get_ylim(info);

% Get positions and heights
v_centers = info.layout.pos(:, 2);
h_centers = info.layout.pos(:, 1);
heights   = info.layout.height;
widths    = info.layout.width;

% Update data of ALL channels, also excluded channels
xmin = min(info.data.time{info.trlop});
xmax = max(info.data.time{info.trlop});
for chanindx = 1:info.nchan
  % Shift the horizontal axis to zero
  hdat = info.data.time{info.trlop} - (xmin + xmax)/2;
  % Scale to length 1 of the new data
  hdat = hdat ./ (xmax - xmin);
  % Scale by availabe plot area
  hdat = hdat .* widths(chanindx);
  % Shift to the horizontal position
  hdat = hdat + h_centers(chanindx);

  % Shift the vertical axis to zero
  vdat = info.data.trial{info.trlop}(chanindx, :); 
  % Scale to length 1 of the new data
  vdat = vdat ./ (ymax - ymin);
  % Scale by availabe plot area
  vdat = vdat .* heights(chanindx);
  % Shift to the vertical position
  vdat = vdat + v_centers(chanindx);
  % Update
  set(info.h_chan(chanindx), {'XData', 'YData'}, {hdat, vdat});
end

% Blank excluded channels
if info.trlsel(info.trlop) % channel is included
  set(info.h_chan(info.chansel), 'Color', 'k');
  set(info.h_chan(~info.chansel), 'Color', 'none');
else % if channel is excluded blank out everything
  set(info.h_chan, 'Color', 'none');
end

% Enable or disable buttons as appropriate
check_buttons_status(info);

% Update title
title(description_trial(info),'interpreter','none');
end

function check_buttons_status(info)
% Updates buttons status if required
%
% Input:
% info - info struct with all the data attched to the GUI

if info.trlop == 1
  info.ui.prev.Enable   = 'off';
  info.ui.prev10.Enable = 'off';
else
  info.ui.prev.Enable   = 'on';
  info.ui.prev10.Enable = 'on';
end
if info.trlop == info.ntrl
  info.ui.next.Enable   = 'off';
  info.ui.next10.Enable = 'off';
else
  info.ui.next.Enable   = 'on';
  info.ui.next10.Enable = 'on';
end
if info.trlop == info.ntrl
  info.ui.excludenext.Enable = 'off';
  info.ui.includenext.Enable = 'off';
else
  info.ui.excludenext.Enable = 'on';
  info.ui.includenext.Enable = 'on';
end
end