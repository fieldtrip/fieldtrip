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
%   cfg.blc       = 'yes';
%   cfg.bpfilter  = [15 25];
%   cfg.bpfreq    =	 'yes';
%   ft_realtime_topography(cfg);

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');
cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% read the header for the first time
hdr = read_header(cfg.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

cfg.channel = channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);

% prepare the layout, also implements channel selection
lay = prepare_layout(cfg);

% determine the size of blocks to process
blocksize = round(cfg.blocksize*hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

% initialize some stuff
cmin = -1;
cmax =  1;
recurz;
% open a new figure
h = figure;

prevSample  = 0;
count       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

  % determine number of samples available in buffer
  hdr = read_header(cfg.headerfile, 'cache', true);

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
    dat = read_data(cfg.datafile, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the power estimation from the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % put the data in a fieldtrip-like raw structure
    data.trial{1} = dat;
    data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
    data.label    = hdr.label(chanindx);
    data.hdr      = hdr;
    data.fsample  = hdr.Fs;

    % apply preprocessing options
    data.trial{1} = preproc(data.trial{1}, data.label, data.fsample, cfg);

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
    tmpcfg            = [];
    tmpcfg.layout     = lay;
    tmpcfg.style      = 'straight';
    tmpcfg.electrodes = 'off';
    tmpcfg.update     = 'no';
    tmpcfg.gridscale  = 35;
    topoplot(tmpcfg, powest);

    c = caxis;
    cmin = min(cmin, c(1));
    cmax = max(cmax, c(2));
    c = [cmin cmax];
    caxis(c);

    drawnow

  end % if enough new samples
end % while true
