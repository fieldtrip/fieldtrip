function artifact_viewer(cfg, artcfg, zval, artval, zindx, inputdata);

% ARTIFACT_VIEWER is a subfunction that reads a segment of data
% (one channel only) and displays it together with the cummulated
% z-value

% Copyright (C) 2004-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_viewer.m,v $
% Revision 1.12  2008/11/18 16:22:42  estmee
% Added cfg.continuous
%
% Revision 1.11  2008/10/07 16:22:12  estmee
% Changed the output of fetch_data and read_data from dat to data.
%
% Revision 1.10  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.9  2006/08/28 08:10:24  jansch
% fixed small bug in number of padding-samples when non-integer, thanks to Jasper
%
% Revision 1.8  2006/01/12 14:20:08  roboos
% new implementation that belongs to artifact_zvalue
%

dat.cfg     = cfg;
dat.artcfg  = artcfg;
if nargin == 5
  % no data is given
  dat.hdr          = read_fcdc_header(cfg.headerfile);
elseif nargin == 6
  % data is given
  dat.hdr          = fetch_header(inputdata); % used name inputdata iso data, because data is already used later in this function
  dat.inputdata    = inputdata; % to be able to get inputdata into h (by guidata)
end
dat.trlop   = 1;
dat.zval    = zval;
dat.artval  = artval;
dat.zindx   = zindx;
dat.stop    = 0;
dat.numtrl  = size(cfg.trl,1);
dat.trialok = zeros(1,dat.numtrl);
for trlop=1:dat.numtrl
  dat.trialok(trlop) = ~any(artval{trlop});
end

h = gcf;
guidata(h,dat);
uicontrol(gcf,'units','pixels','position',[5 5 40 18],'String','stop','Callback',@stop);
uicontrol(gcf,'units','pixels','position',[50 5 25 18],'String','<','Callback',@prevtrial);
uicontrol(gcf,'units','pixels','position',[75 5 25 18],'String','>','Callback',@nexttrial);
uicontrol(gcf,'units','pixels','position',[105 5 25 18],'String','<<','Callback',@prev10trial);
uicontrol(gcf,'units','pixels','position',[130 5 25 18],'String','>>','Callback',@next10trial);
uicontrol(gcf,'units','pixels','position',[160 5 50 18],'String','<artfct','Callback',@prevartfct);
uicontrol(gcf,'units','pixels','position',[210 5 50 18],'String','artfct>','Callback',@nextartfct);

while ishandle(h),
  dat = guidata(h);
  if dat.stop == 0,
    read_and_plot(h);
    uiwait;
  else
    break
  end
end
if ishandle(h)
  close(h);
end

%------------
%subfunctions
%------------

function read_and_plot(h)
vlinecolor = [0 0 0];
dat = guidata(h);
% make a local copy of the relevant variables
trlop = dat.trlop;
zval = dat.zval{trlop};
artval = dat.artval{trlop};
zindx = dat.zindx{trlop};
cfg = dat.cfg;
artcfg = dat.artcfg;
hdr = dat.hdr;
trl = dat.artcfg.trl;
trlpadsmp = round(artcfg.trlpadding*hdr.Fs);
% determine the channel with the highest z-value
[dum, indx] = max(zval);
sgnind = zindx(indx);
iscontinuous = 1;
if isfield(dat, 'inputdata')
  data = fetch_data(dat.inputdata, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
else
  data = read_data(cfg.datafile,   'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
end
% data = preproc(data, channel, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
str = sprintf('trial %3d, channel %s', dat.trlop, hdr.label{sgnind});
fprintf('showing %s\n', str);
% plot z-values in lower subplot
subplot(2,1,2);
cla
hold on
xval = trl(trlop,1):trl(trlop,2);
if trlpadsmp
  sel = 1:trlpadsmp;
  h = plot(xval(sel), zval(sel)); set(h, 'color', [0.5 0.5 1]);  % plot the trialpadding in another color
  sel = trlpadsmp:(length(data)-trlpadsmp);
  plot(xval(sel), zval(sel), 'b');
  sel = (length(data)-trlpadsmp):length(data);
  h = plot(xval(sel), zval(sel)); set(h, 'color', [0.5 0.5 1]);  % plot the trialpadding in another color
  vline(xval(  1)+trlpadsmp, 'color', vlinecolor);
  vline(xval(end)-trlpadsmp, 'color', vlinecolor);
else
  plot(xval, zval, 'b');
end
% draw a line at the threshold level
hline(artcfg.cutoff, 'color', 'r', 'linestyle', ':');
% make the artefact part red
zval(~artval) = nan;
plot(xval, zval, 'r-');
hold off
xlabel('samples');
ylabel('zscore');

% plot data of most aberrant channel in upper subplot
subplot(2,1,1);
cla
hold on
if trlpadsmp
  sel = 1:trlpadsmp;
  h = plot(xval(sel), data(sel)); set(h, 'color', [0.5 0.5 1]);  % plot the trialpadding in another color
  sel = trlpadsmp:(length(data)-trlpadsmp);
  plot(xval(sel), data(sel), 'b');
  sel = (length(data)-trlpadsmp):length(data);
  h = plot(xval(sel), data(sel)); set(h, 'color', [0.5 0.5 1]);  % plot the trialpadding in another color
  vline(xval(  1)+trlpadsmp, 'color', vlinecolor);
  vline(xval(end)-trlpadsmp, 'color', vlinecolor);
else
  plot(xval, data, 'b');
end
data(~artval) = nan;
plot(xval, data, 'r-');
hold off
xlabel('samples');
ylabel('uV or Tesla');
title(str);

function varargout = nexttrial(h, eventdata, handles, varargin)
dat = guidata(h);
if dat.trlop < dat.numtrl,
  dat.trlop = dat.trlop + 1;
end;
guidata(h,dat);
uiresume;

function varargout = next10trial(h, eventdata, handles, varargin)
dat = guidata(h);
if dat.trlop < dat.numtrl - 10,
  dat.trlop = dat.trlop + 10;
else dat.trlop = dat.numtrl;
end;
guidata(h,dat);
uiresume;

function varargout = prevtrial(h, eventdata, handles, varargin)
dat = guidata(h);
if dat.trlop > 1,
  dat.trlop = dat.trlop - 1;
else dat.trlop = 1;
end;
guidata(h,dat);
uiresume;

function varargout = prev10trial(h, eventdata, handles, varargin)
dat = guidata(h);
if dat.trlop > 10,
  dat.trlop = dat.trlop - 10;
else dat.trlop = 1;
end;
guidata(h,dat);
uiresume;

function varargout = nextartfct(h, eventdata, handles, varargin)
dat = guidata(h);
artfctindx = find(dat.trialok == 0);
sel = find(artfctindx > dat.trlop);
if ~isempty(sel)
  dat.trlop = artfctindx(sel(1));
else
  dat.trlop = dat.trlop;
end
guidata(h,dat);
uiresume;

function varargout = prevartfct(h, eventdata, handles, varargin)
dat = guidata(h);
artfctindx = find(dat.trialok == 0);
sel = find(artfctindx < dat.trlop);
if ~isempty(sel)
  dat.trlop = artfctindx(sel(end));
else
  dat.trlop = dat.trlop;
end
guidata(h,dat);
uiresume;

function varargout = stop(h, eventdata, handles, varargin)
dat = guidata(h);
dat.stop = 1;
guidata(h,dat);
uiresume;

function vsquare(x1, x2, c)
abc = axis;
y1 = abc(3);
y2 = abc(4);
x = [x1 x2 x2 x1 x1];
y = [y1 y1 y2 y2 y1];
z = [-1 -1 -1 -1 -1];
h = patch(x, y, z, c);
set(h, 'edgecolor', c)

