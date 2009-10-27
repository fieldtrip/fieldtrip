function rt_headlocalizer(cfg)

% RT_HEADLOCALIZER is an example realtime application for online
% visualization of the head localization coils in a CTF275 system.
%
% Use as
%   rt_headlocalizer(cfg)
% with the following configuration options
%   cfg.template   = string, name of the original data set to be used as template (default = [])
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'last')
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

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isfield(cfg, 'template'),       cfg.template = [];        end
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');

% read the template coil positions
if ~isempty(cfg.template)
  template = read_headshape(cfg.template, 'coordinates', 'dewar');
else
  template = [];
end

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading, only "headloc" type channels are relevant
chanindx = strmatch('headloc', chantype(hdr));

if isempty(chanindx)
  error('the data does not seem to have headlocalization channels');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);

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

  if newsamples>=blocksize

    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample + 1;
      endsample  = prevSample + blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment from buffer
    dat = read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the headlocalization in the CTF system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % put the data in a fieldtrip-like raw structure
    data.trial{1} = dat;
    data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
    data.label    = hdr.label(chanindx);
    data.hdr      = hdr;
    data.fsample  = hdr.Fs;

    % for head coil "n", the labels can be interpreted as
    % HLC00n1 x
    % HLC00n2 y
    % HLC00n3 z
    % HLC00n4 reserved
    % HLC00n5 reserved
    % HLC00n6 reserved
    % HLC00n7 reserved
    % HLC00n8 fit error

    x1i = strmatch('HLC0011', data.label);
    y1i = strmatch('HLC0012', data.label);
    z1i = strmatch('HLC0013', data.label);
    x2i = strmatch('HLC0021', data.label);
    y2i = strmatch('HLC0022', data.label);
    z2i = strmatch('HLC0023', data.label);
    x3i = strmatch('HLC0031', data.label);
    y3i = strmatch('HLC0032', data.label);
    z3i = strmatch('HLC0033', data.label);

    % convert from meter to cm
    coil1 = data.trial{1}([x1i y1i z1i],:) * 100;
    coil2 = data.trial{1}([x2i y2i z2i],:) * 100;
    coil3 = data.trial{1}([x3i y3i z3i],:) * 100;

    figure(1)
    h = get(gca, 'children');
    hold on

    if ~isempty(h)
      % done on every iteration
      delete(h);
    end

    if ~isempty(template)
      % plot the three fiducial positions from the template headcoordinate file
      plot3(template.fid.pnt(1,1), template.fid.pnt(1,2), template.fid.pnt(1,3), 'ko');
      plot3(template.fid.pnt(2,1), template.fid.pnt(2,2), template.fid.pnt(2,3), 'ko');
      plot3(template.fid.pnt(3,1), template.fid.pnt(3,2), template.fid.pnt(3,3), 'ko');
    end

    % plot the coil positons
    plot3(coil1(1,:),coil1(2,:),coil1(3,:), 'r.')
    plot3(coil2(1,:),coil2(2,:),coil2(3,:), 'g.')
    plot3(coil3(1,:),coil3(2,:),coil3(3,:), 'b.')
    % legend({'coil1', 'coil2', 'coil3'}); % this is rather slow according  to the profiler

    str = sprintf('time = %d s\n', round(mean(data.time{1})));
    title(str);
    fprintf(str);

    if isempty(h)
      % only done once
      grid on
      az = 20;
      el = 20;
      xlabel('x (cm)');
      ylabel('y (cm)');
      zlabel('z (cm)');
      set(gca, 'xtick', -30:2.0:30)
      set(gca, 'ytick', -30:2.0:30)
      set(gca, 'ztick', -30:0.5:30) % note the different scaling
      view(az, el);
      % axis square
      axis vis3d
      axis manual
    end

    % force Matlab to update the figure
    drawnow

  end % if enough new samples
end % while true
