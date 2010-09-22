function ft_realtime_signalproxy(cfg)

% FT_REALTIME_SIGNALPROXY creates some random data and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the
% acquisition client to stream data to it. An analysis client can connect
% to read the data upon request. Multiple clients can connect simultaneously,
% each analyzing a specific aspect of the data concurrently.
%
% Use as
%   ft_realtime_signalproxy(cfg)
% with the following configuration options
%   cfg.blocksize            = number, in seconds (default = 1)
%   cfg.channel              = cell-array with channel names
%   cfg.fsample              = sampling frequency
%   cfg.speed                = relative speed at which data is written (default = 1)
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

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'blocksize'),            cfg.blocksize = 1;                                end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = ft_senslabel('eeg1020');            end
if ~isfield(cfg, 'fsample'),              cfg.fsample = 250;                                end % in Hz
if ~isfield(cfg, 'speed'),                cfg.speed = 1 ;                                   end % relative
% set the defaults for filtering
if ~isfield(cfg, 'lpfilter'),             cfg.lpfilter = 'no';                              end
if ~isfield(cfg, 'hpfilter'),             cfg.hpfilter = 'no';                              end
if ~isfield(cfg, 'bpfilter'),             cfg.bpfilter = 'no';                              end
if ~isfield(cfg, 'lpfiltord'),            cfg.lpfiltord = 6;                                end
if ~isfield(cfg, 'hpfiltord'),            cfg.hpfiltord = 6;                                end
if ~isfield(cfg, 'bpfiltord'),            cfg.bpfiltord = 4;                                end
if ~isfield(cfg, 'lpfilttype'),           cfg.lpfilttype = 'but';                           end
if ~isfield(cfg, 'hpfilttype'),           cfg.hpfilttype = 'but';                           end
if ~isfield(cfg, 'bpfilttype'),           cfg.bpfilttype = 'but';                           end
if ~isfield(cfg, 'lpfiltdir'),            cfg.lpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'hpfiltdir'),            cfg.hpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'bpfiltdir'),            cfg.bpfiltdir = 'twopass';                        end
if ~isfield(cfg, 'debug'),                cfg.debug = 'no';                                 end

% translate dataset into datafile+headerfile
cfg.target = checkconfig(cfg.target, 'dataset2files', 'yes');
checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

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

  % wait for a realistic amount of time
  pause(((endsample-begsample+1)/hdr.Fs)/cfg.speed);

  % apply some filters
  if strcmp(cfg.lpfilter, 'yes'),     dat = preproc_lowpassfilter (dat, hdr.Fs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir); end
  if strcmp(cfg.hpfilter, 'yes'),     dat = preproc_highpassfilter(dat, hdr.Fs, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir); end
  if strcmp(cfg.bpfilter, 'yes'),     dat = preproc_bandpassfilter(dat, hdr.Fs, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % from here onward it is specific to writing the data to another stream
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(cfg.debug, 'yes')
    fprintf('sample time = %f, clock time = %f\n', endsample/hdr.Fs, toc(stopwatch));
  end

  if count==1
    % flush the file, write the header and subsequently write the data segment
    write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
  else
    % write the data segment
    write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', true);
  end % if count==1

end % while true
