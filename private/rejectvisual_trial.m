function [chansel, trlsel, cfg] = rejectvisual_trial(cfg, data)

% SUBFUNCTION for ft_rejectvisual

% determine the initial selection of trials
ntrl = length(data.trial);
if isequal(cfg.trials, 'all') % support specification like 'all'
  cfg.trials = 1:ntrl;
end
trlsel  = false(1,ntrl);
trlsel(cfg.trials) = true;

% determine the initial selection of channels
nchan = length(data.label);
cfg.channel = ft_channelselection(cfg.channel, data.label); % support specification like 'all'
chansel = false(1,nchan);
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
for i=1:ntrl
  begsample = nearest(data.time{i}, cfg.latency(1));
  endsample = nearest(data.time{i}, cfg.latency(2));
  data.time{i} = data.time{i}(begsample:endsample);
  data.trial{i} = data.trial{i}(:,begsample:endsample);
end

h = figure;
axis([0 1 0 1]);
axis off

% the info structure will be attached to the figure
% and passed around between the callback functions
if strcmp(cfg.plotlayout,'1col') % hidden config option for plotting trials differently
  info         = [];
  info.ncols   = 1;
  info.nrows   = nchan;
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
  % determine the position of each subplot within the axis
  for row=1:info.nrows
    for col=1:info.ncols
      indx = (row-1)*info.ncols + col;
      if indx>info.nchan
        continue
      end
      info.x(indx)     = (col-0.9)/info.ncols;
      info.y(indx)     = 1 - (row-0.45)/(info.nrows+1);
    end
  end
  info.label = info.data.label;
elseif strcmp(cfg.plotlayout,'square')
  info         = [];
  info.ncols   = ceil(sqrt(nchan));
  info.nrows   = ceil(sqrt(nchan));
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
  % determine the position of each subplot within the axis
  for row=1:info.nrows
    for col=1:info.ncols
      indx = (row-1)*info.ncols + col;
      if indx>info.nchan
        continue
      end
      info.x(indx)     = (col-0.9)/info.ncols;
      info.y(indx)     = 1 - (row-0.45)/(info.nrows+1);
    end
  end
  info.label = info.data.label;
end

info.ui.quit        = uicontrol(h,'units','pixels','position',[  5 5 40 18],'String','quit','Callback',@stop);
info.ui.prev        = uicontrol(h,'units','pixels','position',[ 50 5 25 18],'String','<','Callback',@prev);
info.ui.next        = uicontrol(h,'units','pixels','position',[ 75 5 25 18],'String','>','Callback',@next);
info.ui.prev10      = uicontrol(h,'units','pixels','position',[105 5 25 18],'String','<<','Callback',@prev10);
info.ui.next10      = uicontrol(h,'units','pixels','position',[130 5 25 18],'String','>>','Callback',@next10);
info.ui.bad         = uicontrol(h,'units','pixels','position',[160 5 50 18],'String','bad','Callback',@markbad);
info.ui.good        = uicontrol(h,'units','pixels','position',[210 5 50 18],'String','good','Callback',@markgood);
info.ui.badnext     = uicontrol(h,'units','pixels','position',[270 5 50 18],'String','bad>','Callback',@markbad_next);
info.ui.goodnext    = uicontrol(h,'units','pixels','position',[320 5 50 18],'String','good>','Callback',@markgood_next);
set(gcf, 'WindowButtonUpFcn', @button);
set(gcf, 'KeyPressFcn', @key);

guidata(h,info);

interactive = 1;
while interactive && ishandle(h)
  redraw(h);
  info = guidata(h);
  if info.quit == 0,
    uiwait;
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

function varargout = markbad_next(varargin)
markbad(varargin{:});
next(varargin{:});

function varargout = markgood_next(varargin)
markgood(varargin{:});
next(varargin{:});

function varargout = next(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop < info.ntrl,
  info.trlop = info.trlop + 1;
end;
guidata(h,info);
uiresume;

function varargout = prev(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop > 1,
  info.trlop = info.trlop - 1;
end;
guidata(h,info);
uiresume;

function varargout = next10(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop < info.ntrl - 10,
  info.trlop = info.trlop + 10;
else
  info.trlop = info.ntrl;
end;
guidata(h,info);
uiresume;

function varargout = prev10(h, eventdata, handles, varargin)
info = guidata(h);
info.ltrlop = info.trlop;
if info.trlop > 10,
  info.trlop = info.trlop - 10;
else
  info.trlop = 1;
end;
guidata(h,info);
uiresume;

function varargout = markgood(h, eventdata, handles, varargin)
info = guidata(h);
info.trlsel(info.trlop) = 1;
fprintf(description_trial(info));
title(description_trial(info));
guidata(h,info);
% uiresume;

function varargout = markbad(h, eventdata, handles, varargin)
info = guidata(h);
info.trlsel(info.trlop) = 0;
fprintf(description_trial(info));
title(description_trial(info));
guidata(h,info);
% uiresume;

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
    markgood(h);
  case 'b'
    markbad(h);
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
dx = info.x - x;
dy = info.y - y;
dd = sqrt(dx.^2 + dy.^2);
[d, i] = min(dd);
if d<0.5/max(info.nrows, info.ncols) && i<=info.nchan
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
  str = sprintf('channel %s marked as GOOD\n', info.data.label{info.chanlop});
else
  str = sprintf('channel %s marked as BAD\n', info.data.label{info.chanlop});
end

function str = description_trial(info)
if info.trlsel(info.trlop)
  str = sprintf('trial %d marked as GOOD\n', info.trlop);
else
  str = sprintf('trial %d marked as BAD\n', info.trlop);
end

function redraw(h)
if ~ishandle(h)
  return
end
info = guidata(h);
fprintf(description_trial(info));
cla;
title('');
drawnow
hold on
dat  = info.data.trial{info.trlop};
time = info.data.time{info.trlop};
if ~isempty(info.cfg.alim)
  % use fixed amplitude limits for amplitude scaling
  amax = info.cfg.alim;
else
  % use automatic amplitude limits for scaling, these are different for each trial
  amax = max(max(abs(dat(info.chansel,:))));
end
tmin = time(1);
tmax = time(end);
% scale the time values between 0.1 and 0.9
time = 0.1 + 0.8*(time-tmin)/(tmax-tmin);
% scale the amplitude values between -0.5 and 0.5, offset should not be removed
dat = dat ./ (2*amax);
for row=1:info.nrows
  for col=1:info.ncols
    chanindx = (row-1)*info.ncols + col;
    if chanindx>info.nchan || ~info.chansel(chanindx)
      continue
    end
    % scale the time values for this subplot
    tim = (col-1)/info.ncols + time/info.ncols;
    % scale the amplitude values for this subplot
    amp = dat(chanindx,:)./info.nrows + 1 - row/(info.nrows+1);
    plot(tim, amp, 'k')
  end
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
    set(info.ui.badnext,'Enable', 'off');
    set(info.ui.goodnext,'Enable', 'off');
else
    set(info.ui.badnext,'Enable', 'on');
    set(info.ui.goodnext,'Enable', 'on');
end
text(info.x, info.y, info.label);
title(description_trial(info));
hold off
