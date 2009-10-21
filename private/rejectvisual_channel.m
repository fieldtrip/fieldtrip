function [chansel, trlsel, cfg] = rejectvisual_channel(cfg, data);

% SUBFUNCTION for rejectvisual

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: rejectvisual_channel.m,v $
% Revision 1.4  2007/01/11 13:53:13  roboos
% imlemented cfg.alim, which allows manual specificationh of the amplitude limits in the channel and trial display
%
% Revision 1.3  2007/01/10 11:46:01  roboos
% implemented selection of time window using cfg.latency
%
% Revision 1.2  2006/12/11 10:53:43  roboos
% corrected the internal name of the function (thanks to John Iversen)
%
% Revision 1.1  2006/11/30 13:57:21  roboos
% new implementation, code moved to seperate subfunctions
%

% determine the initial selection of trials and channels
nchan = length(data.label);
ntrl  = length(data.trial);
cfg.channel = channelselection(cfg.channel, data.label);
trlsel  = logical(ones(1,ntrl));
chansel = logical(zeros(1,nchan));
chansel(match_str(data.label, cfg.channel)) = 1;

% compute the offset from the time axes
offset = zeros(ntrl,1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, data.fsample);
end

progress('init', cfg.feedback, 'filtering data');
for i=1:ntrl
  progress(i/ntrl, 'filtering data in trial %d of %d\n', i, ntrl);
  [data.trial{i}, label, time, cfg.preproc] = preproc(data.trial{i}, data.label, data.fsample, cfg.preproc, offset(i));
end
progress('close');

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
info         = [];
info.ncols   = ceil(sqrt(ntrl));
info.nrows   = ceil(sqrt(ntrl));
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
% determine the position of each subplot within the axis
for row=1:info.nrows
  for col=1:info.ncols
    indx = (row-1)*info.ncols + col;
    if indx>info.ntrl
      continue
    end
    info.x(indx)     = (col-0.9)/info.ncols;
    info.y(indx)     = 1 - (row-0.45)/(info.nrows+1);
    info.label{indx} = sprintf('trial %03d', indx);
  end
end

guidata(h,info);

uicontrol(h,'units','pixels','position',[5 5 40 18],'String','quit','Callback',@stop);
uicontrol(h,'units','pixels','position',[50 5 25 18],'String','<','Callback',@prev);
uicontrol(h,'units','pixels','position',[75 5 25 18],'String','>','Callback',@next);
uicontrol(h,'units','pixels','position',[105 5 25 18],'String','<<','Callback',@prev10);
uicontrol(h,'units','pixels','position',[130 5 25 18],'String','>>','Callback',@next10);
uicontrol(h,'units','pixels','position',[160 5 50 18],'String','bad','Callback',@markbad);
uicontrol(h,'units','pixels','position',[210 5 50 18],'String','good','Callback',@markgood);
uicontrol(h,'units','pixels','position',[270 5 50 18],'String','bad>','Callback',@markbad_next);
uicontrol(h,'units','pixels','position',[320 5 50 18],'String','good>','Callback',@markgood_next);
set(gcf, 'WindowButtonUpFcn', @button);
set(gcf, 'KeyPressFcn', @key);

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
if info.chanlop < info.nchan,
  info.chanlop = info.chanlop + 1;
end;
guidata(h,info);
uiresume;

function varargout = prev(h, eventdata, handles, varargin)
info = guidata(h);
if info.chanlop > 1,
  info.chanlop = info.chanlop - 1;
end;
guidata(h,info);
uiresume;

function varargout = next10(h, eventdata, handles, varargin)
info = guidata(h);
if info.chanlop < info.nchan - 10,
  info.chanlop = info.chanlop + 10;
else
  info.chanlop = info.nchan;
end;
guidata(h,info);
uiresume;

function varargout = prev10(h, eventdata, handles, varargin)
info = guidata(h);
if info.chanlop > 10,
  info.chanlop = info.chanlop - 10;
else
  info.chanlop = 1;
end;
guidata(h,info);
uiresume;

function varargout = markgood(h, eventdata, handles, varargin)
info = guidata(h);
info.chansel(info.chanlop) = 1;
fprintf(description_channel(info));
title(description_channel(info));
guidata(h,info);
% uiresume;

function varargout = markbad(h, eventdata, handles, varargin)
info = guidata(h);
info.chansel(info.chanlop) = 0;
fprintf(description_channel(info));
title(description_channel(info));
guidata(h,info);
% uiresume;

function varargout = key(h, eventdata, handles, varargin)
switch lower(eventdata.Key)
  case 'rightarrow'
    next(h);
  case 'leftarrow'
    prev(h);
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
if d<0.5/max(info.nrows, info.ncols) && i<=info.ntrl
  info.trlsel(i) = ~info.trlsel(i); % toggle
  info.trlop = i;
  fprintf(description_trial(info));
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

function str = description_channel(info);
if info.chansel(info.chanlop)
  str = sprintf('channel %s marked as GOOD\n', info.data.label{info.chanlop});
else
  str = sprintf('channel %s marked as BAD\n', info.data.label{info.chanlop});
end

function str = description_trial(info);
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
fprintf(description_channel(info));
cla;
title('');
drawnow
hold on
% determine the maximum value for this channel over all trials
amax = -inf;
tmin =  inf;
tmax = -inf;
for trlindx=find(info.trlsel)
  tmin = min(info.data.time{trlindx}(1)  , tmin);
  tmax = max(info.data.time{trlindx}(end), tmax);
  amax = max(max(abs(info.data.trial{trlindx}(info.chanlop,:))), amax);
end
if ~isempty(info.cfg.alim)
  % use fixed amplitude limits for the amplitude scaling
  amax = info.cfg.alim;
end
for row=1:info.nrows
  for col=1:info.ncols
    trlindx = (row-1)*info.ncols + col;
    if trlindx>info.ntrl || ~info.trlsel(trlindx)
      continue
    end
    % scale the time values between 0.1 and 0.9
    time = info.data.time{trlindx};
    time = 0.1 + 0.8*(time-tmin)/(tmax-tmin);
    % scale the amplitude values between -0.5 and 0.5, offset should not be removed
    dat = info.data.trial{trlindx}(info.chanlop,:);
    dat = dat ./ (2*amax);
    % scale the time values for this subplot
    tim = (col-1)/info.ncols + time/info.ncols;
    % scale the amplitude values for this subplot
    amp = dat./info.nrows + 1 - row/(info.nrows+1);
    plot(tim, amp, 'k')
  end
end
text(info.x, info.y, info.label);
title(description_channel(info));
hold off
