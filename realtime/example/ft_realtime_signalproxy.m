function ft_realtime_signalproxy(cfg)

% FT_REALTIME_SIGNALPROXY creates some random data and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_signalproxy(cfg)
% with the following configuration options
%   cfg.blocksize            = number, in seconds (default = 0.5)
%   cfg.channel              = cell-array with channel names
%   cfg.fsample              = sampling frequency
%   cfg.speed                = relative speed at which data is written (default = 1)
%   cfg.precision            = numeric representation, can be double, single, int32, int16 (default = 'double')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% You can apply some filtering to the random number data to make it
% appear slightly more realistic with
%   cfg.lpfilter      = 'no' or 'yes'  lowpass  filter (default = 'no')
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter (default = 'no')
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter (default = 'no')
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2009-2011, Robert Oostenveld
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
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'blocksize'),            cfg.blocksize = 0.5;                              end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = ft_senslabel('eeg1020');            end
if ~isfield(cfg, 'fsample'),              cfg.fsample = 250;                                end % in Hz
if ~isfield(cfg, 'speed'),                cfg.speed = 1 ;                                   end % relative
% set the defaults for filtering
if ~isfield(cfg, 'onlinefilter'),         cfg.onlinefilter = 'yes';                         end
if ~isfield(cfg, 'lpfilter'),             cfg.lpfilter = 'no';                              end
if ~isfield(cfg, 'hpfilter'),             cfg.hpfilter = 'no';                              end
if ~isfield(cfg, 'bpfilter'),             cfg.bpfilter = 'no';                              end
if ~isfield(cfg, 'lpfiltord'),            cfg.lpfiltord = [];                               end
if ~isfield(cfg, 'hpfiltord'),            cfg.hpfiltord = [];                               end
if ~isfield(cfg, 'bpfiltord'),            cfg.bpfiltord = [];                               end
if ~isfield(cfg, 'lpfilttype'),           cfg.lpfilttype = 'but';                           end
if ~isfield(cfg, 'hpfilttype'),           cfg.hpfilttype = 'but';                           end
if ~isfield(cfg, 'bpfilttype'),           cfg.bpfilttype = 'but';                           end
if ~isfield(cfg, 'lpfiltdir'),            cfg.lpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'hpfiltdir'),            cfg.hpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'bpfiltdir'),            cfg.bpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'debug'),                cfg.debug = 'no';                                 end
if ~isfield(cfg, 'precision'),            cfg.precision = 'double';                         end

% translate dataset into datafile+headerfile
cfg.target = ft_checkconfig(cfg.target, 'dataset2files', 'yes');
ft_checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

hdr = [];
hdr.Fs = cfg.fsample;
hdr.label = cfg.channel;
hdr.nChans = length(cfg.channel);
hdr.nSamples = 0;
hdr.nSamplesPre = 0;
hdr.nTrials = 1;

blocksmp   = round(cfg.blocksize*hdr.Fs);
count      = 0;
prevSample = 0;
stopwatch  = tic;

if strcmp(cfg.onlinefilter, 'yes')
  % Nyquist frequency
  Fn = hdr.Fs/2;
  % oly one type of filter will be applied
  if strcmp(cfg.bpfilter, 'yes')
    Fbp = cfg.bpfreq;
    N   = ft_getopt(cfg, 'bpfiltord', 3);
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);
  elseif strcmp(cfg.lpfilter, 'yes')
    Flp = cfg.lpfreq;
    N   = ft_getopt(cfg, 'lpfiltord', 4);
    [B, A] = butter(N, Flp/Fn);
  elseif strcmp(cfg.hpfilter, 'yes')
    Fhp = cfg.hpfreq;
    N   = ft_getopt(cfg, 'hpfiltord', 4);
    [B, A] = butter(N, Fhp/Fn, 'high');
  end
  % initialize the filter model
  FM = ft_preproc_online_filter_init(B, A, zeros(hdr.nChans, 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true
  
  % increment the number of samples
  hdr.nSamples = hdr.nSamples + blocksmp;
  
  begsample  = prevSample+1;
  endsample  = prevSample+blocksmp;
  
  % remember up to where the data was read
  prevSample  = endsample;
  count       = count + 1;
  fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
  
  % create a random data segment
  dat = randn(hdr.nChans, blocksmp);
  
  switch cfg.precision
    case 'double'
      % keep the data as it is
    case 'single'
      dat = single(dat);
    case 'int32'
      dat = int32(dat);
    case 'int16'
      dat = int16(dat);
    otherwise
      ft_error('unsupported value for cfg.precision');
  end
  
  % wait for a realistic amount of time
  pause(((endsample-begsample+1)/hdr.Fs)/cfg.speed);
  
  % apply some filters
  if strcmp(cfg.onlinefilter, 'yes')
    [FM, dat] = ft_preproc_online_filter_apply(FM, dat);
  else
    if strcmp(cfg.lpfilter, 'yes'),     dat = ft_preproc_lowpassfilter (dat, hdr.Fs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir); end
    if strcmp(cfg.hpfilter, 'yes'),     dat = ft_preproc_highpassfilter(dat, hdr.Fs, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir); end
    if strcmp(cfg.bpfilter, 'yes'),     dat = ft_preproc_bandpassfilter(dat, hdr.Fs, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir); end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to writing the data to another stream
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(cfg.debug, 'yes')
    fprintf('sample time = %f, clock time = %f\n', endsample/hdr.Fs, toc(stopwatch));
  end
  
  if count==1
    % flush the file, write the header and subsequently write the data segment
    ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
  else
    % write the data segment
    ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', true);
  end % if count==1
  
end % while true
