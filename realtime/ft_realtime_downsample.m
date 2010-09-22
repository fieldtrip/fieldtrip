function ft_realtime_downsample(cfg)

% FT_REALTIME_DOWNSAMPLE
%
% Use as
%   ft_realtime_downsample(cfg)
% with the following configuration options
%   cfg.channel              = cell-array, see FT_CHANNELSELECTION (default = 'all')
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
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

% translate dataset into datafile+headerfile
cfg.source = checkconfig(cfg.source, 'dataset2files', 'yes');
cfg.target = checkconfig(cfg.target, 'dataset2files', 'yes');
checkconfig(cfg.source, 'required', {'datafile' 'headerfile'});
checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% read the header for the first time
hdr = read_header(cfg.source.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

targethdr = hdr;
targethdr.Fs = targethdr.Fs/

% define a subset of channels for reading
cfg.channel = channelselection(cfg.channel, hdr.label);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine number of samples available in buffer
  hdr = read_header(cfg.source.headerfile, 'cache', true);

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
    dat = read_data(cfg.source.datafile, 'dataformat', cfg.source.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the downsampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % put the data in a fieldtrip-like raw structure
    data.trial{1} = dat;
    data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
    data.label    = hdr.label(chanindx);
    data.hdr      = hdr;
    data.fsample  = hdr.Fs;

    % apply preprocessing options
    dat = preproc(data.trial{1}, data.label, data.fsample, cfg);
    % do the downsampling
    dat = preproc_resample(dat, hdr.Fs, cfg.fsample, cfg.method);

    if count==1
      % flush the file, write the header and subsequently write the data segment
      write_data(cfg.target.datafile, dat, 'header', targethdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
    else
      % write the data segment
      write_data(cfg.target.datafile, dat, 'header', targethdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
    end

  end % if enough new samples
end % while true

