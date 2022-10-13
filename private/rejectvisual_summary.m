function [chansel, trlsel, cfg] = rejectvisual_summary(cfg, data)

% SUBFUNCTION for ft_rejectvisual

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
% Copyright (C) 2006-2022, Robert Oostenveld, Jan-Mathijs Schoffelen
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

if ft_platform_supports('matlabversion', '2014b', inf)
  % make the figure large enough to hold stuff and place it centered at the top of the
  % screen where it would appear by default
  pos = get(groot, 'DefaultFigurePosition'); % groot was introduced in R2014b
  w = 800;
  h = 600;
  x = pos(1) + pos(3)/2 - w/2;
  y = pos(2) + pos(4)/2 - h/2;
  cfg.position = ft_getopt(cfg, 'position', [x y w h]);
end

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

% % select the specified latency window from the data
% % here it is done BEFORE filtering and metric computation
% for i=1:ntrl
%   begsample = nearest(data.time{i}, cfg.latency(1));
%   endsample = nearest(data.time{i}, cfg.latency(2));
%   data.time{i} = data.time{i}(begsample:endsample);
%   data.trial{i} = data.trial{i}(:, begsample:endsample);
% end
if ischar(cfg.latency)
  cfg.latency(1) = min(cellfun(@min, data.time));
  cfg.latency(2) = max(cellfun(@max, data.time));
end

% compute the offset from the time axes
offset = zeros(ntrl, 1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, fsample);
end

% set up guidata info
info                = [];
info.data           = data;
info.cfg            = cfg;
info.metric         = cfg.metric;
info.previousmetric = 'none';
info.level          = nan(nchan, ntrl);
info.ntrl           = ntrl;
info.nchan          = nchan;
info.trlsel         = trlsel;
info.chansel        = chansel;
info.fsample        = fsample;
info.offset         = offset;
info.quit           = 0;

if isfield(cfg, 'neighbours')
  % prepare the neighbours, load from disk, make channel selection, etc.
  tmpcfg = keepfields(cfg, {'neighbours'});
  info.cfg.neighbours = ft_prepare_neighbours(tmpcfg, data);
  % creates a NxN Boolean matrix that describes whether channels are connected as neighbours
  info.cfg.connectivity = channelconnectivity(info.cfg, data);
else
  info.cfg.connectivity = [];
end

h = create_figure(info);

% Compute initial metric...
compute_metric(h);

while ishandle(h)
  redraw(h);
  info = guidata(h);
  if info.quit == 0
    uiwait(h);
  else
    chansel = info.chansel;
    trlsel  = info.trlsel;
    cfg     = info.cfg;
    delete(h);
    break;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = create_figure(info)
% Creates the GUI and plots the data for the first time.
% All additional changes are done by manipulating the properties of the
% specific graphical objects
h = open_figure(keepfields(info.cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));

% define three axes
info.axes(1) = axes('position', [0.100 0.650 0.375 0.300]);  % summary
info.axes(2) = axes('position', [0.575 0.650 0.375 0.300]);  % channels
info.axes(3) = axes('position', [0.100 0.250 0.375 0.300]);  % trials
% callback function (for toggling trials/channels) is set later, so that
% the user cannot try to toggle trials while nothing is present on the
% plots

% set up radio buttons for choosing metric
bgcolor = get(h, 'color');

metriclist = {'var' 'std' 'min' 'max' 'maxabs' 'range' 'kurtosis' '1/var' 'zvalue' 'maxzvalue' 'neighbexpvar' 'neighbcorr' 'neighbstdratio' 'db'};

uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.50 0.24 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Metric used:'); % text string for metric
uicontrol(h, 'Units', 'normalized', 'position', [0.725 0.51 0.20 0.04], 'Style', 'popupmenu', 'backgroundcolor', bgcolor, 'string', metriclist,       'HandleVisibility', 'off', 'callback', @change_metric); % popup menu for metric
uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.46 0.24 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Plot trial:');
info.plottrltxt = uicontrol(h, 'Units', 'normalized', 'position', [0.725 0.47 0.20 0.04], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @display_trial); % editbox for trial plotting

% editboxes for manually specifying which channels/trials to toggle on or off
uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.42 0.24 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Toggle trial:');
uicontrol(h, 'Units', 'normalized', 'position', [0.725 0.43 0.20 0.04], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @toggle_trials, 'tag', 'edit_toggle_trials');
uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.38 0.24 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Toggle channel:');
uicontrol(h, 'Units', 'normalized', 'position', [0.725 0.39 0.20 0.04], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @toggle_channels, 'tag', 'edit_toggle_channels');

info.excludetrllbl  = uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.34 0.24 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', sprintf('Rejected trials: %i/%i', sum(info.trlsel==0), info.ntrl));
info.excludetrltxt  = uicontrol(h, 'Units', 'normalized', 'position', [0.605 0.30 0.45 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'Enable', 'Inactive');
info.excludechanlbl = uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.26 0.28 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', sprintf('Rejected channels: %i/%i', sum(info.chansel==0), info.nchan));
info.excludechantxt = uicontrol(h, 'Units', 'normalized', 'position', [0.605 0.22 0.45 0.04], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'Enable', 'Inactive');

% instructions
instructions = sprintf('Drag the mouse over the channels or trials you wish to exclude');
uicontrol(h, 'Units', 'normalized', 'position', [0.575 0.125 0.35 0.08], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', instructions, 'FontWeight', 'bold', 'ForegroundColor', 'r');

% logbox
info.output_box = uicontrol(h, 'Units', 'normalized', 'position', [0.00 0.00 1.00 0.1], 'Style', 'edit', 'HorizontalAlignment', 'left', 'Max', 3, 'Min', 1, 'Enable', 'inactive', 'FontName', get(0, 'FixedWidthFontName'), 'FontSize', 9, 'ForegroundColor', [0 0 0], 'BackgroundColor', [1 1 1]);

% quit button
uicontrol(h, 'Units', 'normalized', 'position', [0.825 0.11 0.10 0.04], 'string', 'quit', 'callback', @quit);

% store the updated info in the figure
guidata(h, info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compute_metric(h)
info = guidata(h);
tmp = info.level(info.chansel, info.trlsel);
if isequal(info.metric, info.previousmetric) && all(~isnan(tmp(:)))
  % there is no reason to recompute the metric
  return
end

if contains(info.metric, 'zvalue')
  % cellmean and cellstd (see FT_DENOISE_PCA) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(info.nchan, 1);
  runss  = zeros(info.nchan, 1);
  runnum = 0;
  for i=1:info.ntrl
    % not entirely sure whether info.data.time{i} is correct, so I am constructing it here on the fly
    dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc);
    runsum = runsum + nansum(dat, 2);
    runss  = runss  + nansum(dat.^2, 2);
    runnum = runnum + sum(isfinite(dat), 2);
  end
  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
else
  mval = [];
  sd   = [];
end

update_log(info.output_box, sprintf('Computing %s...', info.metric));
ft_progress('init', info.cfg.feedback, sprintf('computing %s', info.metric));
level = nan(info.nchan, info.ntrl);
for i=1:info.ntrl
  ft_progress(i/info.ntrl, 'computing %s %d of %d\n', i, info.metric, info.ntrl);
  % not entirely sure whether info.data.time{i} is correct, so I am constructing it here on the fly
  dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc);
  val = artifact_level(dat, info.metric, mval, sd, info.cfg.connectivity);
  if isvector(val)
    level(:,i) = val;
  else
    level(:,i) = nanmedian(val);
  end
end

ft_progress('close');
update_log(info.output_box, 'Done.');

info.level = level;
info.previousmetric = info.metric;
guidata(h, info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%
function redraw(h)
info  = guidata(h);
% work with a copy of the data
level = info.level;

[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, info.chansel, info.trlsel, info.metric);

% make the three figures
if gcf~=h, figure(h); end
% datacursormode on;

set(h, 'CurrentAxes', info.axes(1))
cla(info.axes(1));
switch info.cfg.viewmode
  case {'remove', 'toggle'}
    tmp = level;
    tmp(~info.chansel, :) = nan;
    tmp(:, ~info.trlsel)  = nan;
    imagesc(tmp, 'AlphaData', ~isnan(tmp));
    caxis([min(tmp(:)) max(tmp(:))]);
  case 'hide'
    imagesc(level(info.chansel==1, info.trlsel==1));
    if ~all(info.trlsel)
      set(info.axes(1), 'Xtick', []);
    end
    if ~all(info.chansel)
      set(info.axes(1), 'Ytick', []);
    end
end % switch
axis ij;
% colorbar;
title(info.cfg.method);
ylabel('channel number');
xlabel('trial number');

set(h, 'CurrentAxes', info.axes(2))
cla(info.axes(2));
switch info.cfg.viewmode
  case 'remove'
    plot(maxperchan(info.chansel==1),     find(info.chansel==1), '.');
    xmin = min(maxperchan);
    xmax = max(maxperchan);
    ymin = 1;
    ymax = info.nchan;
  case 'toggle'
    plot(maxperchan_all(info.chansel==1), find(info.chansel==1), '.');
    hold on;
    plot(maxperchan_all(info.chansel==0), find(info.chansel==0), 'o');
    hold off;
    xmin = min(maxperchan_all);
    xmax = max(maxperchan_all);
    ymin = 1;
    ymax = info.nchan;
  case 'hide'
    xmin = min(maxperchan);
    xmax = max(maxperchan);
    ymin = 1;
    ymax = sum(info.chansel==1);
    plot(maxperchan(info.chansel==1), 1:ymax, '.');
    if ~all(info.chansel)
      set(info.axes(2), 'Ytick', []);
    end
end % switch
% don't try to rescale the axes if they are empty
if any(info.chansel) && any(info.trlsel)
  % ensure that the horizontal and vertical range increase, also when negative
  % see https://github.com/fieldtrip/fieldtrip/issues/1150
  axis(fixrange([0.8*xmin 1.2*xmax ymin-0.5 ymax+0.5]));
end
axis ij;
set(info.axes(2), 'ButtonDownFcn', @toggle_visual);  % needs to be here; call to axis resets this property
ylabel('channel number');

set(h, 'CurrentAxes', info.axes(3))
cla(info.axes(3));
switch info.cfg.viewmode
  case 'remove'
    plot(find(info.trlsel==1), maxpertrl(info.trlsel==1), '.');
    xmin = 1;
    xmax = info.ntrl;
    ymin = min(maxpertrl);
    ymax = max(maxpertrl);
  case 'toggle'
    plot(find(info.trlsel==1), maxpertrl_all(info.trlsel==1), '.');
    hold on;
    plot(find(info.trlsel==0), maxpertrl_all(info.trlsel==0), 'o');
    hold off;
    xmin = 1;
    xmax = info.ntrl;
    ymin = min(maxpertrl_all);
    ymax = max(maxpertrl_all);
  case 'hide'
    xmin = 1;
    xmax = sum(info.trlsel==1);
    ymin = min(maxpertrl);
    ymax = max(maxpertrl);
    plot(1:xmax, maxpertrl(info.trlsel==1), '.');
    if ~all(info.trlsel)
      set(info.axes(3), 'Xtick', []);
    end
end % switch
% don't try to rescale the axes if they are empty
if any(info.chansel) && any(info.trlsel)
  % ensure that the horizontal and vertical range increase, also when negative
  % see https://github.com/fieldtrip/fieldtrip/issues/1150
  axis(fixrange([xmin-0.5 xmax+0.5 0.8*ymin 1.2*ymax]));
end
set(info.axes(3), 'ButtonDownFcn', @toggle_visual);  % needs to be here; call to axis resets this property
xlabel('trial number');

% put excluded trials/channels in their respective edit boxes
set(info.excludechanlbl, 'string', sprintf('Channels to exclude: %i/%i', sum(info.chansel==0), info.nchan));
set(info.excludetrllbl, 'string', sprintf('Trials to exclude: %i/%i', sum(info.trlsel==0), info.ntrl));

if ~all(info.trlsel)
  excludetrltxt = sprintf('%d, ', find(~info.trlsel));
  excludetrltxt = excludetrltxt(1:end-2);
  set(info.excludetrltxt, 'String', excludetrltxt, 'FontAngle', 'normal', 'ButtonDownFcn', @toggle_trials_fill);
else
  set(info.excludetrltxt, 'String', 'No trials to exclude', 'FontAngle', 'italic');
end

if ~all(info.chansel)
  if false % isfield(info.data, 'label')
    excludechantxt = sprintf('%s, ', info.data.label{~info.chansel});
  else
    excludechantxt = sprintf('%d, ', find(~info.chansel));
  end
  excludechantxt = excludechantxt(1:end-2);
  set(info.excludechantxt, 'String', excludechantxt, 'FontAngle', 'normal', 'ButtonDownFcn', @toggle_channels_fill);
else
  set(info.excludechantxt, 'String', 'No channels to exclude', 'FontAngle', 'italic');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_trials(h, eventdata)
info = guidata(h);
% process input from the "toggle trials" textbox
rawtrls = get(h, 'string');
set(h, 'string', '');
if ~isempty(rawtrls)
  spltrls = regexp(rawtrls, '\s+', 'split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
else
  update_log(info.output_box, sprintf('Please enter one or more trials'));
  uiresume;
  return;
end
try
  toggle = trls;
  info.trlsel(toggle) = ~info.trlsel(toggle);
catch
  update_log(info.output_box, sprintf('ERROR: Trial value too large!'));
  uiresume;
  return;
end
% recalculate the metric
compute_metric(h)
guidata(h, info);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_trials_fill(h,eventdata)
info = guidata(h);
excludetrltxt = sprintf('%d, ', find(~info.trlsel));
excludetrltxt = excludetrltxt(1:end-2);
set(findobj('tag', 'edit_toggle_trials'), 'string', excludetrltxt)
set(info.plottrltxt, 'string', excludetrltxt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_channels(h, eventdata)
% process input from the "toggle channels" textbox
info = guidata(h);
rawchans = get(h, 'string');
set(h, 'string', '');
if ~isempty(rawchans)
  splchans = regexp(rawchans, '\s+', 'split');
  chans = zeros(1, length(splchans));
  % determine whether identifying channels via number or label
  [junk, junk, junk, procchans] = regexp(rawchans, '([A-Za-z]+|[0-9]{4, })');
  clear junk;
  if isempty(procchans)
    % if using channel numbers
    for n = 1:length(splchans)
      chans(n) = str2num(splchans{n});
    end
  else
    % if using channel labels
    for n = 1:length(splchans)
      try
        chans(n) = find(ismember(info.data.label, splchans(n)));
      catch
        update_log(info.output_box, sprintf('ERROR: Please ensure the channel name is correct (case-sensitive)!'));
        uiresume;
        return;
      end
    end
  end
else
  update_log(info.output_box, sprintf('Please enter one or more channels'));
  uiresume;
  return;
end
try
  toggle = chans;
  info.chansel(toggle) = ~info.chansel(toggle);
catch
  update_log(info.output_box, sprintf('ERROR: Channel value too large!'));
  uiresume;
  return
end
% recalculate the metric
compute_metric(h)
guidata(h, info);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_channels_fill(h, eventdata)
info = guidata(h);
excludechantxt = sprintf('%d, ', find(~info.chansel));
excludechantxt = excludechantxt(1:end-2);
set(findobj('tag', 'edit_toggle_channels'), 'string', excludechantxt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_visual(h, eventdata)
% copied from FT_SELECT_BOX, but without the waitforbuttonpress command since here it  is triggered by the ButtonDown event
point1 = get(gca, 'CurrentPoint');    % button down detected
finalRect = rbbox;                    % return figure units
point2 = get(gca, 'CurrentPoint');    % button up detected
point1 = point1(1, 1:2);              % extract x and y
point2 = point2(1, 1:2);
x = sort([point1(1) point2(1)]);
y = sort([point1(2) point2(2)]);

g     = get(gca, 'Parent');
info  = guidata(g);

[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(info.level, info.chansel, info.trlsel, info.metric);

switch gca
  case info.axes(1)
    % visual selection in the summary plot is not supported
    
  case info.axes(2)
    % the visual selection was made in the channels plot
    switch info.cfg.viewmode
      case 'toggle'
        chanlabels = 1:info.nchan;
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan_all(chanlabels)' >= x(1) & ...
          maxperchan_all(chanlabels)' <= x(2);
        info.chansel(toggle) = ~info.chansel(toggle);
        
      case 'remove'
        chanlabels     = 1:info.nchan;
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan(chanlabels)' >= x(1) & ...
          maxperchan(chanlabels)' <= x(2);
        info.chansel(toggle) = false;
        
      case 'hide'
        chanlabels = 1:sum(info.chansel==1);
        [junk, origchanlabels] = find(info.chansel==1);
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan(origchanlabels)' >= x(1) & ...
          maxperchan(origchanlabels)' <= x(2);
        info.chansel(origchanlabels(toggle)) = false;
    end
    
  case info.axes(3)
    % the visual selection was made in the trials plot
    switch info.cfg.viewmode
      case 'toggle'
        trllabels = 1:info.ntrl;
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl_all(trllabels) >= y(1) & ...
          maxpertrl_all(trllabels) <= y(2);
        info.trlsel(toggle) = ~info.trlsel(toggle);
        
      case 'remove'
        trllabels = 1:info.ntrl;
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl(trllabels) >= y(1) & ...
          maxpertrl(trllabels) <= y(2);
        info.trlsel(toggle) = false;
        
      case 'hide'
        trllabels = 1:sum(info.trlsel==1);
        [junk, origtrllabels] = find(info.trlsel==1);
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl(origtrllabels) >= y(1) & ...
          maxpertrl(origtrllabels) <= y(2);
        info.trlsel(origtrllabels(toggle)) = false;
    end
    
    
end % switch gca

% recalculate the metric
compute_metric(h);
guidata(h, info);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quit(h, eventdata)
info = guidata(h);
info.quit = 1;
guidata(h, info);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function change_metric(h, eventdata)
info = guidata(h);
%info.metric = get(eventdata.NewValue, 'string');
info.metric = h.String{h.Value};
guidata(h, info);
compute_metric(h);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_rejected(h, eventdata)
info = guidata(h);
toggle = get(h, 'value');
if toggle == 0
  info.cfg.viewmode = 'remove';
else
  info.cfg.viewmode = 'toggle';
end
guidata(h, info);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_log(h, new_text)
new_text        = [datestr(now, 13) '# ' new_text];
curr_text       = get(h, 'string');
size_curr_text  = size(curr_text, 2);
size_new_text   = size(new_text, 2);
if size_curr_text > size_new_text
  new_text        = [new_text blanks(size_curr_text-size_new_text)];
else
  curr_text   = [curr_text repmat(blanks(size_new_text-size_curr_text), size(curr_text, 1), 1)];
end
set(h, 'String', [new_text; curr_text]);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, chansel, trlsel, metric)
if ismember(metric, {'var', 'std', 'max', 'maxabs', 'range', 'kurtosis', '1/var', 'zvalue', 'maxzvalue', 'neighbstdratio'})
  % take the maximum
  extreme = @max;
elseif ismember(metric, {'min', 'neighbexpvar', 'neighbcorr'})
  % take the minimum
  extreme = @min;
end
% determine the extreme value
maxperchan_all = extreme(level, [], 2);
maxpertrl_all  = extreme(level, [], 1);
% determine the extreme value over the remaining selection
level(~chansel, :) = nan;
level(:, ~trlsel)  = nan;
maxperchan = extreme(level, [], 2);
maxpertrl  = extreme(level, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_trial(h, eventdata)
info = guidata(h);
update_log(info.output_box, 'Making multiplot of individual trials ...');
rawtrls = get(h, 'string');
set(h, 'string', '');
if isempty(rawtrls)
  return;
else
  spltrls = regexp(rawtrls, '\s+', 'split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
end
cfg_mp = [];
% disable hashing of input data (speeds up things)
cfg_mp.trackcallinfo = 'no';
cfg_mp.layout  = info.cfg.layout;
cfg_mp.channel = info.data.label(info.chansel);
cfg_mp.dataname = info.cfg.dataname;
cfg_mp.ylim     = info.cfg.ylim;
currfig = gcf;
for n = 1:length(trls)
  % ft_multiplotER should be able to make the selection, but fails due to http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2978
  % that bug is hard to fix, hence it is solved here with a work-around
  cfg_sd = [];
  cfg_sd.trials = trls(n);
  tmpdata = ft_selectdata(cfg_sd, info.data);
  
  figure()
  cfg_mp.interactive = 'yes';
  ft_multiplotER(cfg_mp, tmpdata);
  title(sprintf('Trial %i', trls(n)));
end
figure(currfig);
update_log(info.output_box, 'Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range = fixrange(range)
% ensure that the horizontal and vertical range always increase
% also when both are negative, or both are the same, or both are or zero
% see https://github.com/fieldtrip/fieldtrip/issues/1150
if range(1)>range(2)
  % swap them
  range([1 2]) = range([2 1]);
elseif range(1)==range(2)
  if range(1)==0
    % move them a little bit apart from each other
    range([1 2]) = range([1 2]) + [-eps +eps];
  else
    % move them a little bit apart from each other, scaled with the value of interest
    range([1 2]) = range([1 2]) + [-eps +eps]*range(1);
  end
end
if range(3)>range(4)
  % swap them
  range([3 4]) = range([4 3]);
elseif range(3)==range(4)
  if range(3)==0
    % move them a little bit apart from each other
    range([3 4]) = range([3 4]) + [-eps +eps];
  else
    % move them a little bit apart from each other, scaled with the value of interest
    range([3 4]) = range([3 4]) + [-eps +eps]*range(3);
  end
end
