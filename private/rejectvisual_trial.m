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

if ischar(cfg.latency)
  cfg.latency(1) = min(cellfun(@min, data.time));
  cfg.latency(2) = max(cellfun(@max, data.time));
end


% the info structure will be attached to the figure
% and passed around between the callback functions
info         = [];
info.chanlop = 1;
info.trlop   = 1;
info.ltrlop  = 0;
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

h = figure;
axis off
axis([minx maxx miny maxy]);

info.ui.quit        = uicontrol(h,'units','pixels','position',[  5 5 40 18],'String','quit','Callback',@stop);
info.ui.prev        = uicontrol(h,'units','pixels','position',[ 50 5 25 18],'String','<','Callback',@prev);
info.ui.next        = uicontrol(h,'units','pixels','position',[ 75 5 25 18],'String','>','Callback',@next);
info.ui.prev10      = uicontrol(h,'units','pixels','position',[105 5 25 18],'String','<<','Callback',@prev10);
info.ui.next10      = uicontrol(h,'units','pixels','position',[130 5 25 18],'String','>>','Callback',@next10);
info.ui.exclude     = uicontrol(h,'units','pixels','position',[160 5 70 18],'String','exclude','Callback',@markexclude);
info.ui.include     = uicontrol(h,'units','pixels','position',[230 5 70 18],'String','include','Callback',@markinclude);
info.ui.excludenext = uicontrol(h,'units','pixels','position',[310 5 70 18],'String','exclude >','Callback',@markexclude_next);
info.ui.includenext = uicontrol(h,'units','pixels','position',[380 5 70 18],'String','include >','Callback',@markinclude_next);
set(gcf, 'WindowButtonUpFcn', @button);
set(gcf, 'KeyPressFcn', @key);

guidata(h,info);

interactive = 1;
while interactive && ishandle(h)
  redraw(h);
  info = guidata(h);
  if info.quit == 0
    uiwait(h);
  else
    chansel = info.chansel;
    trlsel = info.trlsel;
    delete(h);
    break
  end
end % while interactive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = markexclude_next(varargin)
markexclude(varargin{:});
next(varargin{:});

function varargout = markinclude_next(varargin)
markinclude(varargin{:});
next(varargin{:});

function varargout = next(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop < info.ntrl
  info.trlop = info.trlop + 1;
end
guidata(h,info);
uiresume;

function varargout = prev(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop > 1
  info.trlop = info.trlop - 1;
end
guidata(h,info);
uiresume;

function varargout = next10(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop < info.ntrl - 10
  info.trlop = info.trlop + 10;
else
  info.trlop = info.ntrl;
end
guidata(h,info);
uiresume;

function varargout = prev10(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop > 10
  info.trlop = info.trlop - 10;
else
  info.trlop = 1;
end
guidata(h,info);
uiresume;

function varargout = markinclude(h, eventdata, handles, varargin)
info = guidata(h);
info.trlsel(info.trlop) = 1;
fprintf(description_trial(info));
title(description_trial(info));
guidata(h,info);
uiresume;

function varargout = markexclude(h, eventdata, handles, varargin)
info = guidata(h);
info.trlsel(info.trlop) = 0;
fprintf(description_trial(info));
title(description_trial(info));
guidata(h,info);
uiresume;

function varargout = key(h, eventdata, handles, varargin)
info = guidata(h);
switch lower(eventdata.Key)
  case 'rightarrow'
    if info.trlop ~= info.ntrl
      next(h);
    else
      fprintf('at last trial\n');
    end
  case 'leftarrow'
    if info.trlop ~= 1
      prev(h);
    else
      fprintf('at first trial\n');
    end
  case 'g'
    markinclude(h);
  case 'b'
    markexclude(h);
  case 'q'
    stop(h);
  otherwise
    fprintf('unknown key pressed\n');
end

function varargout = button(h, eventdata, handles, varargin)
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
  info.chansel(i) = ~info.chansel(i); % toggle
  info.chanlop = i;
  fprintf(description_channel(info));
  guidata(h,info);
  uiresume;
else
  fprintf('button clicked\n');
  return
end

function varargout = stop(h, eventdata, handles, varargin)
info = guidata(h);
info.quit = 1;
guidata(h,info);
uiresume;

function str = description_channel(info)
if info.chansel(info.chanlop)
  str = sprintf('channel %s marked to INCLUDE\n', info.data.label{info.chanlop});
else
  str = sprintf('channel %s marked to EXCLUDE\n', info.data.label{info.chanlop});
end

function str = description_trial(info)
if info.trlsel(info.trlop)
  str = sprintf('trial %d marked to INCLUDE\n', info.trlop);
else
  str = sprintf('trial %d marked to EXCLUDE\n', info.trlop);
end

function redraw(h)
if ~ishandle(h)
  return
end
info = guidata(h);
cla;
title('');
drawnow
hold on

if ~isnumeric(info.cfg.ylim)
  ymin = [];
  ymax = [];
  for chanindx=1:info.nchan
    if ~info.chansel(chanindx)
      % do not consider excluded channels in vertical scaling estimate
      continue
    end
    dat = info.data.trial{info.trlop}(chanindx,:);
    ymin = min(ymin, min(dat));
    ymax = max(ymax, max(dat));
  end
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

for chanindx=1:info.nchan
  tim = info.data.time{info.trlop};
  dat = info.data.trial{info.trlop}(chanindx,:);
  if info.chansel(chanindx) && info.trlsel(info.trlop)
    color = 'k';
  else
    color = 'none';
  end
  ft_plot_vector(tim, dat, 'hpos', info.layout.pos(chanindx,1), 'vpos', info.layout.pos(chanindx,2), 'width', info.layout.width(chanindx), 'height', info.layout.height(chanindx), 'vlim', [ymin ymax], 'box', istrue(info.cfg.box), 'color', color, 'label', info.layout.label{chanindx}, 'labelpos', 'lowerleft');
end

% enable or disable buttons as appropriate
if info.trlop == 1
  set(info.ui.prev,   'Enable', 'off');
  set(info.ui.prev10, 'Enable', 'off');
else
  set(info.ui.prev,   'Enable', 'on');
  set(info.ui.prev10, 'Enable', 'on');
end
if info.trlop == info.ntrl
  set(info.ui.next,   'Enable', 'off');
  set(info.ui.next10, 'Enable', 'off');
else
  set(info.ui.next,   'Enable', 'on');
  set(info.ui.next10, 'Enable', 'on');
end
if info.ltrlop == info.trlop && info.trlop == info.ntrl
  set(info.ui.excludenext,'Enable', 'off');
  set(info.ui.includenext,'Enable', 'off');
else
  set(info.ui.excludenext,'Enable', 'on');
  set(info.ui.includenext,'Enable', 'on');
end
title(description_trial(info));
hold off
