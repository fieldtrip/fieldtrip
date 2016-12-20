function ft_realtime_topography(cfg)

% FT_REALTIME_TOPOGRAPHY reads continuous data from a file or from a data stream,
% estimates the power and plots the scalp topography in real time.
%
% Use as
%   ft_realtime_topography(cfg)
% with the following configuration options
%   cfg.blocksize            = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.overlap              = number, amojunt of overlap between chunks (default = 0 seconds)
%   cfg.layout               = specification of the layout, see FT_PREPARE_LAYOUT
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
% To stop this realtime function, you have to press Ctrl-C
%
% Example use
%   cfg           = [];
%   cfg.dataset   = 'PW02_ingnie_20061212_01.ds';
%   cfg.layout    = 'CTF151.lay';
%   cfg.channel   = 'MEG';
%   cfg.blocksize = 0.5;
%   cfg.overlap   = 0.25;
%   cfg.demean    = 'yes';
%   cfg.bpfilter  = [15 25];
%   cfg.bpfreq    =	 'yes';
%   ft_realtime_topography(cfg);

% Copyright (C) 2008, Robert Oostenveld
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

% set the defaults
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% read the header for the first time
hdr = ft_read_header(cfg.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);

% prepare the layout, also implements channel selection
lay = ft_prepare_layout(cfg);

% determine the size of blocks to process
blocksize = round(cfg.blocksize*hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

% initialize some stuff
cmin = -1;
cmax =  1;
clear recurz
recurz; % initialize the persistent variables

% open a new figure
h = figure;

prevSample  = 0;
count       = 0;

lay = ft_prepare_layout(cfg);
[laysel, datsel] = match_str(lay.label, hdr.label);
% get the 2D position of the channels
x = lay.pos(laysel,1);
y = lay.pos(laysel,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'cache', true);

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples>=(blocksize-overlap)

    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample + 1;
      endsample  = prevSample + blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end

    % this allows overlapping data segments
    if overlap && (begsample>overlap)
      begsample = begsample - overlap;
      endsample = endsample - overlap;
    end

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment
    dat = ft_read_data(cfg.datafile, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the power estimation from the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % put the data in a fieldtrip-like raw structure
    data          = [];
    data.trial{1} = dat;
    data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
    data.label    = hdr.label(chanindx);
    data.hdr      = hdr;
    data.fsample  = hdr.Fs;

    % apply preprocessing options
    data = ft_preprocessing(cfg, data);

    % estimate power
    powest = sum(data.trial{1}.^2, 2);

    if ~ishandle(h)
      % re-initialize some stuff
      cmin = -1;
      cmax =  1;
      % open a new figure
      h = figure;
    end

    % compute z-transformed
    powest = recurz(powest);

    % plot the topography
    ft_plot_topo(x, y, powest(datsel), 'outline', lay.outline, 'mask', lay.mask);
    hold on
    plot(x, y, 'k.');
    hold off

    c = caxis;
    cmin = min(cmin, c(1));
    cmax = max(cmax, c(2));
    c = [cmin cmax];
    caxis(c);

    drawnow

  end % if enough new samples
end % while true


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time = offset2time(offset, fsample, nsamples)
time = (offset + (0:(nsamples-1)))/fsample;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION recursive computation of z-transformed data by means of persistent variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = recurz(x)

persistent n
persistent s
persistent ss

if nargin==0 || isempty(x)
  % re-initialize
  n  = [];
  s  = [];
  ss = [];
  return
end

if isempty(n)
  n = 1;
else
  n = n + 1;
end

if isempty(s)
  s = x;
else
  s = s + x;
end

if isempty(ss)
  ss = x.^2;
else
  ss = ss + x.^2;
end

if n==1
  % standard deviation cannot be computed yet
  z = zeros(size(x));
elseif all(s(:)==ss(:))
  % standard deviation is zero anyway
  z = zeros(size(x));
else
  % compute standard deviation and z-transform of the input data
  sd = sqrt((ss - (s.^2)./n) ./ (n-1));
  z  = (x-s/n)./ sd;
end
