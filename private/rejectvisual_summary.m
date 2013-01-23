function [chansel, trlsel, cfg] = rejectvisual_summary(cfg, data)

% REJECTVISUAL_SUMMARY:  subfunction for ft_rejectvisual

% determine the initial selection of trials and channels
nchan = length(data.label);
ntrl  = length(data.trial);
cfg.channel = ft_channelselection(cfg.channel, data.label);
trlsel  = true(1,ntrl);
chansel = false(1,nchan);
chansel(match_str(data.label, cfg.channel)) = 1;

% compute the sampling frequency from the first two timepoints
fsample = 1/mean(diff(data.time{1}));


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
  offset(i) = time2offset(data.time{i}, fsample);
end


% set up guidata info
h = figure();
info = [];
info.data       = data;
info.cfg        = cfg;
info.metric     = cfg.metric;
info.ntrl       = ntrl;
info.nchan      = nchan;
info.trlsel     = trlsel;
info.chansel    = chansel;
info.fsample    = fsample;
info.offset     = offset;
info.quit       = 0;
guidata(h,info);

% set up display
interactive = 1;
info = guidata(h);

% make the figure large enough to hold stuff
set(h,'Position',[50 350 800 500]);

% define three plots
info.axes(1) = axes('position',[0.100 0.650 0.375 0.300]);  % summary plot
info.axes(2) = axes('position',[0.575 0.650 0.375 0.300]);  % channels
info.axes(3) = axes('position',[0.100 0.250 0.375 0.300]);  % trials
% callback function (for toggling trials/channels) is set later, so that
% the user cannot try to toggle trials while nothing is present on the
% plots

% set up radio buttons for choosing metric
bgcolor = get(h,'color');
g = uibuttongroup('Position',[0.525 0.275 0.375 0.250 ],'bordertype','none','backgroundcolor',bgcolor);
r(1) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 7/7 0.40 0.15 ],'Style','radio','backgroundcolor',bgcolor,'string','var','HandleVisibility','off');
r(2) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 6/7 0.40 0.15 ],'Style','radio','backgroundcolor',bgcolor,'String','min','HandleVisibility','off');
r(3) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 5/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','max','HandleVisibility','off');
r(4) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 4/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','maxabs','HandleVisibility','off');
r(5) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 3/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','range','HandleVisibility','off');
r(6) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 2/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','kurtosis','HandleVisibility','off');
r(7) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 1/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','1/var','HandleVisibility','off');
r(8) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 0/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','zvalue','HandleVisibility','off');
r(9) = uicontrol('Units','normalized','parent',g,'position',[ 0.0 -1/7 0.40 0.15 ],'Style','Radio','backgroundcolor',bgcolor,'String','maxzvalue','HandleVisibility','off');
% pre-select appropriate metric, if defined
set(g,'SelectionChangeFcn',@change_metric);
for i=1:length(r)
    if strcmp(get(r(i),'string'), cfg.metric)
        set(g,'SelectedObject',r(i));
    end
end

% editboxes for manually specifying which channels/trials to turn off
uicontrol(h,'Units','normalized','position',[0.64 0.44 0.14 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string','Toggle trial #:');
uicontrol(h,'Units','normalized','position',[0.64 0.40 0.12 0.05],'Style','edit','HorizontalAlignment','left','backgroundcolor',[1 1 1],'callback',@toggle_trials);
uicontrol(h,'Units','normalized','position',[0.64 0.31 0.14 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string','Toggle channel:');
uicontrol(h,'Units','normalized','position',[0.64 0.27 0.12 0.05],'Style','edit','HorizontalAlignment','left','backgroundcolor',[1 1 1],'callback',@toggle_channels);

% editbox for trial plotting
                  uicontrol(h,'Units','normalized','position',[0.500 0.165 0.10 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string','Plot trial:');
info.plottrltxt = uicontrol(h,'Units','normalized','position',[0.580 0.170 0.12 0.05],'Style','edit','HorizontalAlignment','left','backgroundcolor',[1 1 1],'callback',@display_trial);

info.badtrllbl  = uicontrol(h,'Units','normalized','position',[0.795 0.44 0.195 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string',sprintf('Rejected trials: %i/%i',sum(info.trlsel==0),info.ntrl));
info.badtrltxt  = uicontrol(h,'Units','normalized','position',[0.795 0.3975 0.23 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'));
info.badchanlbl = uicontrol(h,'Units','normalized','position',[0.795 0.31 0.195 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string',sprintf('Rejected channels: %i/%i',sum(info.chansel==0),info.nchan));
info.badchantxt = uicontrol(h,'Units','normalized','position',[0.795 0.2625 0.23 0.05],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'));

% instructions
instructions = sprintf('Drag the mouse over the channels/trials you wish to reject.');
uicontrol(h,'Units','normalized','position',[0.625 0.50 0.35 0.065],'Style','text','HorizontalAlignment','left','backgroundcolor',get(h,'color'),'string',instructions,'FontWeight','bold','ForegroundColor','b');

% "show rejected" button
% ui_tog = uicontrol(h,'Units','normalized','position',[0.55 0.200 0.25 0.05],'Style','checkbox','backgroundcolor',get(h,'color'),'string','Show rejected?','callback',@toggle_rejected);
% if strcmp(cfg.viewmode, 'toggle')
%     set(ui_tog,'value',1);
% end

% logbox
info.output_box = uicontrol(h,'Units','normalized','position',[0.00 0.00 1.00 0.15],'Style','edit','HorizontalAlignment','left','Max',3,'Min',1,'Enable','inactive','FontName',get(0,'FixedWidthFontName'),'FontSize',9,'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]);

% quit button
uicontrol(h,'Units','normalized','position',[0.80 0.175 0.10 0.05],'string','quit','callback',@quit);

guidata(h, info);

% disable trial plotting if cfg.layout not present
if ~isfield(info.cfg,'layout')
  set(info.plottrltxt,'Enable','off');
  update_log(info.output_box,sprintf('NOTE: "cfg.layout" parameter required for trial plotting!'));
end

% Compute initial metric...
compute_metric(h);

while interactive && ishandle(h)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_metric(h)
info = guidata(h);
update_log(info.output_box,'Computing metric...');
ft_progress('init', info.cfg.feedback, 'computing metric');
level = zeros(info.nchan, info.ntrl);
if strcmp(info.metric,'zvalue') || strcmp(info.metric, 'maxzvalue')
    % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, 
    % but they were too memory-intensive
    runsum=zeros(info.nchan, 1);
    runss=zeros(info.nchan,1);
    runnum=0;
    for i=1:info.ntrl
        [dat] = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i},2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
        dat(info.chansel==0,:) = nan;
        runsum=runsum+sum(dat,2);
        runss=runss+sum(dat.^2,2);
        runnum=runnum+size(dat,2);
    end
    mval=runsum/runnum;
    sd=sqrt(runss/runnum - (runsum./runnum).^2);
end
for i=1:info.ntrl
  ft_progress(i/info.ntrl, 'computing metric %d of %d\n', i, info.ntrl);
  [dat, label, time, info.cfg.preproc] = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i},2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
  dat(info.chansel==0,:) = nan;
  switch info.metric
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
    case 'zvalue'
      level(:,i) = mean( ( dat-repmat(mval,1,size(dat,2)) )./repmat(sd,1,size(dat,2)) ,2);
    case 'maxzvalue'
      level(:,i) = max( ( dat-repmat(mval,1,size(dat,2)) )./repmat(sd,1,size(dat,2)) , [], 2);
    otherwise
      error('unsupported method');
  end
end
ft_progress('close');
update_log(info.output_box,'Done.');
% if metric calculated with channels off, then turning that channel on
% won't show anything. Because of this, keep track of which channels were
% off when calculated, so know when to recalculate.
info.metric_chansel = info.chansel;
% % reinsert the data for the selected channels
% dum = nan(info.nchan, info.ntrl);
% dum(info.chansel,:) = level;
% origlevel = dum;
% clear dum
% % store original levels for later
% info.origlevel = origlevel;
info.level = level;
guidata(h,info);

function redraw(h)
info  = guidata(h);
% work with a copy of the data
level = info.level;
level(~info.chansel,:) = nan;
level(:,~info.trlsel)  = nan;

[maxperchan maxpertrl maxperchan_all maxpertrl_all] = set_maxper(level, info.chansel, info.trlsel);

% make the three figures
if gcf~=h, figure(h); end
%datacursormode on;

set(h,'CurrentAxes',info.axes(1))
cla(info.axes(1));
% imagesc(level(chansel, trlsel));
imagesc(level);
axis xy;
% colorbar;
title(info.cfg.method);
ylabel('channel number');
xlabel('trial number');

set(h,'CurrentAxes',info.axes(2))
cla(info.axes(2));
set(info.axes(2),'ButtonDownFcn',@toggle_visual);
plot(maxperchan(info.chansel==1),     find(info.chansel==1), '.');
if strcmp(info.cfg.viewmode, 'toggle') && (sum(info.chansel==0) > 0)
  hold on;
  plot(maxperchan_all(info.chansel==0), find(info.chansel==0), 'o');
  hold off;
end
abc = axis;
axis([abc(1:2) 0 info.nchan]); % have to use 0 as lower limit because ylim([1 1]) (i.e. the single-channel case) is invalid
set(info.axes(2),'ButtonDownFcn',@toggle_visual);  % needs to be here; call to axis resets this property
ylabel('channel number');

set(h,'CurrentAxes',info.axes(3))
cla(info.axes(3));
plot(find(info.trlsel==1), maxpertrl(info.trlsel==1), '.');
if strcmp(info.cfg.viewmode, 'toggle') && (sum(info.trlsel==0) > 0)
  hold on;
  plot(find(info.trlsel==0), maxpertrl_all(info.trlsel==0), 'o');
  hold off;
end
abc = axis; axis([1 info.ntrl abc(3:4)]);
set(info.axes(3),'ButtonDownFcn',@toggle_visual);  % needs to be here; call to axis resets this property
xlabel('trial number');

% put rejected trials/channels in their respective edit boxes
set(info.badchanlbl,'string',sprintf('Rejected channels: %i/%i',sum(info.chansel==0),info.nchan));
set(info.badtrllbl, 'string',sprintf('Rejected trials: %i/%i',sum(info.trlsel==0), info.ntrl));
if ~isempty(find(info.trlsel==0, 1))
  set(info.badtrltxt,'String',num2str(find(info.trlsel==0)),'FontAngle','normal');
else
  set(info.badtrltxt,'String','No trials rejected','FontAngle','italic');
end
if ~isempty(find(info.chansel==0, 1))
  if isfield(info.data,'label')
    chanlabels = info.data.label(info.chansel==0);
    badchantxt = '';
      for i=find(info.chansel==0)
        if ~isempty(badchantxt)
          badchantxt = [badchantxt ', ' info.data.label{i} '(' num2str(i) ')'];
        else
          badchantxt = [info.data.label{i} '(' num2str(i) ')'];
        end
      end
    set(info.badchantxt,'String',badchantxt,'FontAngle','normal');
  else
    set(info.badtrltxt,'String',num2str(find(info.chansel==0)),'FontAngle','normal');
  end
else
  set(info.badchantxt,'String','No channels rejected','FontAngle','italic');
end

function toggle_trials(h, eventdata)
info = guidata(h);
% extract trials from string
rawtrls = get(h,'string');
if ~isempty(rawtrls)
  spltrls = regexp(rawtrls,'\s+','split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
else
  update_log(info.output_box,sprintf('Please enter one or more trials'));
  uiresume;
  return;
end
toggle = trls;
try
  info.trlsel(toggle) = ~info.trlsel(toggle);
catch
  update_log(info.output_box,sprintf('ERROR: Trial value too large!'));
end

guidata(h, info);
uiresume;

% process input from the "toggle channels" textbox
function toggle_channels(h, eventdata)
info = guidata(h);
rawchans = get(h,'string');
if ~isempty(rawchans)
  splchans = regexp(rawchans,'\s+','split');
  chans = zeros(1,length(splchans));
  % determine whether identifying channels via number or label
  [junk junk junk procchans] = regexp(rawchans,'([A-Za-z]+|[0-9]{4,})');
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
        chans(n) = find(ismember(info.data.label,splchans(n)));
      catch
        update_log(info.output_box,sprintf('ERROR: Please ensure the channel name is correct (case-sensitive)!'));
        uiresume;
        return;
      end
    end
  end    
else
  update_log(info.output_box,sprintf('Please enter one or more channels'));
  uiresume;
  return;
end
toggle = chans;
try
  info.chansel(toggle) = ~info.chansel(toggle);
  % if levels data from channel being toggled was calculated from another
  % metric, recalculate the metric
  if info.metric_chansel(toggle) == 0
    compute_metric(h)
  end
catch
  update_log(info.output_box,sprintf('ERROR: Channel value too large!'));
end
guidata(h, info);

uiresume;

function toggle_visual(h, eventdata)
% copied from select2d, without waitforbuttonpress command
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
x = sort([point1(1) point2(1)]);
y = sort([point1(2) point2(2)]);

g     = get(gca,'Parent');
info  = guidata(g);

[maxperchan maxpertrl] = set_maxper(info.level, info.chansel, info.trlsel);

% toggle channels
if gca == info.axes(2)
%   if strcmp(info.cfg.viewmode, 'toggle')
    chanlabels = 1:info.nchan;
%   else
%     perchan    = max(info.level,[],2);
%     maxperchan = perchan(info.chansel==1);
%     chanlabels = find(info.chansel==1);
%   end
  toggle = find( ...
    chanlabels >= y(1) & ...
    chanlabels <= y(2) & ...
    maxperchan(:)' >= x(1) & ...
    maxperchan(:)' <= x(2));
  info.chansel(toggle) = 0;
  % if levels data from channel being toggled was calculated from another
  % metric, recalculate the metric
  if info.metric_chansel(toggle) == 0
      compute_metric(h)
  end
% toggle trials
elseif gca == info.axes(3)
%   if strcmp(info.cfg.viewmode, 'toggle')
    trllabels = 1:info.ntrl;
%   else
%     pertrl    = max(info.level,[],1);
%     maxpertrl = pertrl(info.trlsel==1);
%     trllabels = find(info.trlsel==1);
%   end
  toggle = find( ...
    trllabels >= x(1) & ...
    trllabels <= x(2) & ...
    maxpertrl(:)' >= y(1) & ...
    maxpertrl(:)' <= y(2));
  info.trlsel(toggle) = 0;
end
guidata(h, info);
uiresume;

% function display_trial(h, eventdata)
% info = guidata(h);
% rawtrls = get(h,'string');
% if ~isempty(rawtrls)
%   spltrls = regexp(rawtrls,' ','split');
%   trls = [];
%   for n = 1:length(spltrls)
%     trls(n) = str2num(cell2mat(spltrls(n)));
%   end
% else
%   update_log(info.output_box,sprintf('Please enter one or more trials'));
%   uiresume;
%   return;
% end
% if all(trls==0)
%   % use visual selection
%   update_log(info.output_box,sprintf('make visual selection of trials to be plotted seperately...'));
%   [x, y] = select2d;
%   maxpertrl  = max(info.origlevel,[],1);
%   toggle = find(1:ntrl>=x(1) & ...
%     1:ntrl<=x(2) & ...
%     maxpertrl(:)'>=y(1) & ...
%     maxpertrl(:)'<=y(2));
% else
%   toggle = trls;
% end
% for i=1:length(trls)
%   figure
%   % the data being displayed here is NOT filtered
%   %plot(data.time{toggle(i)}, data.trial{toggle(i)}(chansel,:));
%   tmp = info.data.trial{toggle(i)}(info.chansel,:);
%   tmp = tmp - repmat(mean(tmp,2), [1 size(tmp,2)]);
%   plot(info.data.time{toggle(i)}, tmp);
%   title(sprintf('trial %d', toggle(i)));
% end
  
function quit(h, eventdata)
info = guidata(h);
info.quit = 1;
guidata(h, info);
uiresume;

function change_metric(h, eventdata)
info = guidata(h);
info.metric = get(eventdata.NewValue, 'string');
guidata(h, info);
compute_metric(h);
uiresume;

function toggle_rejected(h, eventdata)
info = guidata(h);
toggle = get(h,'value');
if toggle == 0
  info.cfg.viewmode = 'remove';
else
  info.cfg.viewmode = 'toggle';
end
guidata(h, info);
uiresume;

function update_log(h, new_text)
new_text        = [datestr(now,13) '# ' new_text];
curr_text       = get(h, 'string');
size_curr_text  = size(curr_text,2);
size_new_text   = size(new_text,2);
if size_curr_text > size_new_text
  new_text        = [new_text blanks(size_curr_text-size_new_text)];
else
  curr_text   = [curr_text repmat(blanks(size_new_text-size_curr_text), size(curr_text,1), 1)];
end
set(h, 'String', [new_text; curr_text]);
drawnow;

function [maxperchan,maxpertrl,varargout] = set_maxper(level,chansel,trlsel)
level_all = level;
% determine the maximum value
level(~chansel,:) = nan;
level(:,~trlsel)  = nan;

maxperchan = max(level,[],2);
maxpertrl  = max(level,[],1);
maxperchan_all = max(level_all,[],2);
maxpertrl_all  = max(level_all,[],1);
varargout(1) = {maxperchan_all};
varargout(2) = {maxpertrl_all};

function display_trial(h, eventdata)
info = guidata(h);
rawtrls = get(h,'string');
if isempty(rawtrls)
  return;
else
  spltrls = regexp(rawtrls,'\s+','split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
end
cfg_mp = [];
cfg_mp.layout  = info.cfg.layout;
cfg_mp.channel = info.data.label(info.chansel);
currfig = gcf;
for n = 1:length(trls)
  figure()
  cfg_mp.trials = trls(n);
  cfg_mp.interactive = 'yes';
  ft_multiplotER(cfg_mp, info.data);
  title(sprintf('Trial %i',trls(n)));
end
figure(currfig);
return;
