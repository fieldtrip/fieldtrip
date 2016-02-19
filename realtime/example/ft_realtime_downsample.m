function ft_realtime_downsample(cfg)

% FT_REALTIME_DOWNSAMPLE reads realtime data from one buffer and writes it after downsampling
% to another buffer.
%
% Use as
%   ft_realtime_downsample(cfg)
% with the following configuration options
%   cfg.channel              = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.decimation           = integer, downsampling factor (default = 1, no downsampling)
%   cfg.order                = interger, order of butterworth lowpass filter (default = 4)
%   cfg.cutoff               = double, cutoff frequency of lowpass filter (default = 0.8*Nyquist-freq.)
%
% The source of the data is configured as
%   cfg.source.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.source.datafile      = string
%   cfg.source.headerfile    = string
%   cfg.source.eventfile     = string
%   cfg.source.dataformat    = string, default is determined automatic
%   cfg.source.headerformat  = string, default is determined automatic
%   cfg.source.eventformat   = string, default is determined automatic
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2008, Robert Oostenveld
% Copyright (C) 2010, Stefan Klanke
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
if ~isfield(cfg, 'source'),               cfg.source = [];                                  end
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.source, 'headerformat'),  cfg.source.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.source, 'dataformat'),    cfg.source.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'minblocksize'),         cfg.minblocksize = 0;                             end % in seconds
if ~isfield(cfg, 'maxblocksize'),         cfg.maxblocksize = 1;                             end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = 'all';                              end

if ~isfield(cfg, 'decimation'),           cfg.decimation = 1;                               end
if ~isfield(cfg, 'order'),                cfg.order = 4;                                    end

if cfg.decimation < 1 || cfg.decimation~=round(cfg.decimation)
  error 'Decimation factor must be integer quantity >= 1';
end

% translate dataset into datafile+headerfile
cfg.source = ft_checkconfig(cfg.source, 'dataset2files', 'yes');
cfg.target = ft_checkconfig(cfg.target, 'dataset2files', 'yes');
ft_checkconfig(cfg.source, 'required', {'datafile' 'headerfile'});
ft_checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% read the header for the first time
hdr = ft_read_header(cfg.source.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

targethdr = hdr;
targethdr.Fs = targethdr.Fs/cfg.decimation;

% Calculate cutoff frequency expressed in Nyquist freq of original signal
if ~isfield(cfg, 'cutoff')
  cutoff = 0.8/cfg.decimation;
else
  cutoff = cfg.cutoff / (0.5*hdr.Fs);
  if cutoff >= 1
    error('Cutoff frequency too high');
  end
end
[B,A] = butter(cfg.order, cutoff, 'low');

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

minblocksmp = round(cfg.minblocksize*hdr.Fs);
minblocksmp = max(minblocksmp, 1);
maxblocksmp = round(cfg.maxblocksize*hdr.Fs);
prevSample  = 0;
count       = 0;

% create filter model and downsampling model
FM = online_filter_init(B, A, zeros(nchan, 1));
DM = online_downsample_init(cfg.decimation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.source.headerfile, 'cache', true);

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples>=minblocksmp

    % determine the samples to copy from source to target
    begsample = prevSample+1;
    endsample = prevSample + min(newsamples, maxblocksmp);

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment
    dat = ft_read_data(cfg.source.datafile, 'dataformat', cfg.source.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the downsampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [FM, dat] = online_filter_apply(FM, dat);
    [DM, xd] = online_downsample_apply(DM, dat);

    if count==1
      % flush the file, write the header and subsequently write the data segment
      ft_write_data(cfg.target.datafile, dat, 'header', targethdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
    else
      % write the data segment
      ft_write_data(cfg.target.datafile, dat, 'header', targethdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
    end

  end % if enough new samples
end % while true

