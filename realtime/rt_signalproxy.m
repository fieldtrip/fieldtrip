function rt_signalproxy(cfg)

% RT_SIGNALPROXY creates some random data and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the
% acquisition client to stream data to it. An analysis client can connect
% to read the data upon request. Multiple clients can connect simultaneously,
% each analyzing a specific aspect of the data concurrently.
%
% Use as
%   rt_signalproxy(cfg)
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
% $Log: rt_signalproxy.m,v $
% Revision 1.6  2009/05/01 15:41:21  roboos
% consistent handling of cfg.speed
%
% Revision 1.5  2009/05/01 14:08:30  roboos
% updated documentation
%
% Revision 1.4  2009/03/30 11:59:45  roboos
% fixed typo
%
% Revision 1.3  2009/03/30 11:59:26  roboos
% fixed documentation, added elapsed and signal time to the fprintf statement
%
% Revision 1.2  2009/01/21 21:00:08  roboos
% added more elaborate filter options
% added cfg.speed to run slower/faster than realtime
%
% Revision 1.1  2009/01/20 15:43:24  roboos
% new function meant for testing
%

% set the defaults
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'blocksize'),            cfg.blocksize = 1;                                end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = senslabel('eeg1020');               end
if ~isfield(cfg, 'fsample'),              cfg.fsample = 250;                                end
if ~isfield(cfg, 'bpfreq'),               cfg.bpfreq = [];                                  end
if ~isfield(cfg, 'speed'),                cfg.speed = 1 ;                                   end
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
t0         = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true

  % increment the number of samples
  hdr.nSamples = hdr.nSamples + blocksmp;

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples>=blocksmp

    begsample  = prevSample+1;
    endsample  = prevSample+blocksmp;

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d, stime = %f, etime = %f\n', count, begsample, endsample, hdr.nSamples/hdr.Fs, etime(clock, t0));

    % create a random data segment
    dat = randn(hdr.nChans, blocksmp);

    % apply some filters
    if strcmp(cfg.lpfilter, 'yes'),     dat = preproc_lowpassfilter (dat, hdr.Fs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir); end
    if strcmp(cfg.hpfilter, 'yes'),     dat = preproc_highpassfilter(dat, hdr.Fs, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir); end
    if strcmp(cfg.bpfilter, 'yes'),     dat = preproc_bandpassfilter(dat, hdr.Fs, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to writing the data to another stream
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if count==1
      % flush the file, write the header and subsequently write the data segment
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
    else
      % write the data segment
      write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', true);
    end % if count==1

    % wait for a realistic amount of time
    pause(((endsample-begsample+1)/hdr.Fs)/cfg.speed);

  end % if enough new samples
end % while true
