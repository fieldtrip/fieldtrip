function ft_realtime_oddball(cfg)

% FT_REALTIME_ODDBALL is an realtime application that computes an online
% average for a standard and deviant condition. The ERPs/ERFs are plotted,
% together with the difference as t-values. It should work both for EEG and
% MEG, as long as there are two triggers present
%
% Use as
%   ft_realtime_oddball(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.trialfun   = string with the trial function
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

% Copyright (C) 2008-2012, Robert Oostenveld
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

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end % jump to end of file at initialization

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% these are used by the GUI callbacks
clear global chansel chanindx vaxis hdr
global chansel chanindx vaxis hdr

b1clicked = false;
b2clicked = false;
chansel = 1; % this is the subselection out of chanindx
vaxis = [
  -6 6
  -3 3
  ];

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample = 0;
end
count = 0;

f1 = nan;

% initialize the timelock cell-array, each cell will hold the average in one condition
timelock = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
  
  % determine latest header and event information
  event     = ft_read_event(cfg.dataset, 'minsample', prevSample+1);  % only consider events that are later than the data processed sofar
  hdr       = ft_read_header(cfg.dataset, 'cache', true);             % the trialfun might want to use this, but it is not required
  cfg.event = event;                                               % store it in the configuration, so that it can be passed on to the trialfun
  cfg.hdr   = hdr;                                                 % store it in the configuration, so that it can be passed on to the trialfun
  
  % evaluate the trialfun, note that the trialfun should not re-read the events and header
  fprintf('evaluating ''%s'' based on %d events\n', cfg.trialfun, length(event));
  trl = feval(cfg.trialfun, cfg);
  
  % the code below assumes that the 4th column of the trl matrix contains the condition index
  % set the default condition to one if no condition index was given
  if size(trl,1)>0 && size(trl,2)<4
    trl(:,4) = 1;
  end
  
  fprintf('processing %d trials\n', size(trl,1));
  
  for trllop=1:size(trl,1)
    
    begsample = trl(trllop,1);
    endsample = trl(trllop,2);
    offset    = trl(trllop,3);
    condition = trl(trllop,4);  % it is important that the 4th column is returned with the condition number
    
    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d, condition = %d\n', count, begsample, endsample, condition);

    while (hdr.nSamples*hdr.nTrials < endsample)
      % wait until all data up to the endsample has arrived
      hdr = ft_read_header(cfg.headerfile, 'cache', true);
    end
    
    % read the selected data segment from the buffer
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the processing of the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % apply some preprocessing options
    dat = ft_preproc_lowpassfilter(dat, hdr.Fs, 45);
    dat = ft_preproc_baselinecorrect(dat, 1, -offset);
    
    % put the data in a fieldtrip-like raw structure
    data.trial{1} = dat;
    data.time{1}  = offset2time(offset, hdr.Fs, endsample-begsample+1);
    data.label    = hdr.label(chanindx);
    data.hdr      = hdr;
    data.fsample  = hdr.Fs;
    
    if length(timelock)<condition || isempty(timelock{condition})
      % this is the first occurence of this condition, initialize an empty timelock structure
      timelock{condition}.label   = data.label;
      timelock{condition}.time    = data.time{1};
      timelock{condition}.avg     = [];
      timelock{condition}.var     = [];
      timelock{condition}.dimord  = 'chan_time';
      nchans   = size(data.trial{1}, 1);
      nsamples = size(data.trial{1}, 2);
      % the following elements are for the cumulative computation
      timelock{condition}.n   = 0;                          % number of trials
      timelock{condition}.s   = zeros(nchans, nsamples);    % sum
      timelock{condition}.ss  = zeros(nchans, nsamples);    % sum of squares
    end
    
    % add the new data to the accumulated data
    timelock{condition}.n  = timelock{condition}.n  + 1;
    timelock{condition}.s  = timelock{condition}.s  + data.trial{1};
    timelock{condition}.ss = timelock{condition}.ss + data.trial{1}.^2;
    
    % compute the average and variance on the fly
    timelock{condition}.avg = timelock{condition}.s ./ timelock{condition}.n;
    timelock{condition}.var = (timelock{condition}.ss - (timelock{condition}.s.^2)./timelock{condition}.n) ./ (timelock{condition}.n-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward the GUI is constructed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
      
      if ~ishandle(f1)
        close all;
        f1 = figure;
        clear u1 u2
        clear p1 p2
        clear c1
        set(f1, 'resizeFcn', 'clear u1 u2 p1 p2 c1 b1 b2')
      end
      
      if ~exist('p1')
        p1 = subplot(2,1,1);
      end
      
      if ~exist('p2')
        p2 = subplot(2,1,2);
      end
      
      if ~exist('c1')
        pos = [0.75 0.93 0.1 0.05];
        c1 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_channel, 'BackgroundColor', 'white');
        set(c1, 'position', pos);
        set(c1, 'string', chanindx(chansel));
        set(c1, 'tag', 'c1');
      end
      
      if ~exist('u1')
        pos = get(p1, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u1 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u1, 'position', pos);
        set(u1, 'string', num2str(vaxis(1,2)));
        set(u1, 'tag', 'u1');
      end
      
      if ~exist('u2')
        pos = get(p2, 'position'); % link the position to the subplot
        pos(1) = pos(1)-0.1;
        pos(2) = pos(2)-0.05;
        pos(3) = 0.1;
        pos(4) = 0.05;
        u2 = uicontrol('style', 'edit', 'units', 'normalized', 'callback', @update_axis, 'BackgroundColor', 'white');
        set(u2, 'position', pos);
        set(u2, 'string', num2str(vaxis(2,2)));
        set(u2, 'tag', 'u1');
      end
      
      if ~exist('b1')
        pos = [0.75 0.01 0.1 0.05];
        b1 = uicontrol('style', 'pushbutton', 'units', 'normalized', 'callback', 'evalin(''caller'', ''b1clicked = true'')');
        set(b1, 'position', pos);
        set(b1, 'string', 'reset');
        set(b1, 'tag', 'b1');
      end
      
      if ~exist('b2')
        pos = [0.88 0.01 0.1 0.05];
        b2 = uicontrol('style', 'pushbutton', 'units', 'normalized', 'callback', 'evalin(''caller'', ''b2clicked = true'')');
        set(b2, 'position', pos);
        set(b2, 'string', 'quit');
        set(b2, 'tag', 'b2');
      end
      
    end % try
    
    if b1clicked
      timelock = {};
      try, cla(p1); end
      try, cla(p2); end
      b1clicked = false;
    end
    
    if b2clicked
      return
      b2clicked = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward the data is plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
      
      if length(timelock)>1
        sel = ~cellfun(@isempty, timelock);
        sel = find(sel, 2, 'first');
        if length(sel)~=2
          break
        end
        standard  = timelock{sel(1)};
        deviant   = timelock{sel(2)};
        tscore    = (deviant.avg - standard.avg) ./ sqrt(standard.var./standard.n + deviant.var./deviant.n);
        time      = standard.time; % deviant is the same
        
        if exist('p1')
          subplot(p1)
          cla
          hold on
          hs = plot(time, standard.avg(chansel,:), 'b-');
          hd = plot(time, deviant.avg(chansel,:), 'r-');
          set(hs, 'lineWidth', 1.5)
          set(hd, 'lineWidth', 1.5)
          grid on
          axis([time(1) time(end) vaxis(1,1) vaxis(1,2)])
          legend(sprintf('standard (n=%d)', standard.n), sprintf('deviant (n=%d)', deviant.n));
          xlabel('time (s)');
          ylabel('amplitude (uV)');
          title(sprintf('channel "%s"', hdr.label{chanindx(chansel)}));
        end
        
        if exist('p2')
          subplot(p2)
          cla
          hold on
          ht = plot(time, tscore(chansel,:), 'g-');
          set(ht, 'lineWidth', 1.5)
          grid on
          axis([time(1) time(end) vaxis(2,1) vaxis(2,2)])
          legend('difference');
          xlabel('time (s)');
          ylabel('t-score (a.u.)');
        end
        
      end % two conditions are available
      
    end % try
    
    % force matlab to redraw the figure
    drawnow
    
  end % looping over new trials
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)
offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_channel(h, varargin)
global chansel chanindx hdr
val = abs(str2num(get(h, 'string')));
val = max(1, min(val, length(chanindx)));
if ~isempty(val)
  switch get(h, 'tag')
    case 'c1'
      chansel = val;
      set(h, 'string', num2str(val));
      fprintf('switching to channel "%s"', hdr.label{chanindx(chansel)});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_axis(h, varargin)
global vaxis
val = abs(str2num(get(h, 'string')));
if ~isempty(val)
  switch get(h, 'tag')
    case 'u1'
      vaxis(1,:) = [-val val];
    case 'u2'
      vaxis(2,:) = [-val val];
  end
end
