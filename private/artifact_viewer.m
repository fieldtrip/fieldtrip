function artifact_viewer(cfg, artcfg, zval, artval, zindx, inputdata)

% ARTIFACT_VIEWER is a subfunction that reads a segment of data
% (one channel only) and displays it together with the cumulated
% z-value

% Copyright (C) 2004-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
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

if ishandle(cfg)
  % get the input variables from the handle's guidata
  hx     = cfg; clear cfg;
  tmp    = guidata(hx);
  cfg    = tmp.cfg;
  cfg.trl = tmp.trl;
  artcfg = tmp.artcfg;
  zval   = tmp.zval;
  artval = tmp.artval;
  zindx  = tmp.zindx;
  if ~isempty(tmp.data)
    inputdata = tmp.data;
    hdr       = ft_fetch_header(inputdata);
  else
    hdr       = ft_read_header(cfg.headerfile);
  end
  % make a new figure
  figure;
else
  if nargin == 5
    hdr = ft_read_header(cfg.headerfile);
  elseif nargin == 6
    hdr = ft_fetch_header(inputdata); % used name inputdata iso data, because data is already used later in this function
  end
end

dat.cfg     = cfg;
dat.artcfg  = artcfg;
if exist('inputdata', 'var') dat.inputdata = inputdata; end
dat.trlop   = 1;
dat.zval    = zval;
dat.artval  = artval;
dat.zindx   = zindx;
dat.stop    = 0;
dat.numtrl  = size(cfg.trl,1);
dat.trialok = zeros(1,dat.numtrl);
dat.hdr     = hdr;
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
  data = ft_fetch_data(dat.inputdata, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
else
  data = ft_read_data(cfg.datafile,   'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
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

