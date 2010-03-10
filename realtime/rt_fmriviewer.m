function rt_fmriviewer(cfg)

% RT_FMRIVIEWER
%
% Use as
%   rt_fmriviewer(cfg)
% with the following configuration options
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

if ~isfield(cfg, 'headerfile'),    	cfg.headerfile = 'buffer://localhost:1972'; end
if ~isfield(cfg, 'datafile'),      	cfg.datafile = cfg.headerfile; end
if ~isfield(cfg, 'eventfile'),    	cfg.eventfile = cfg.headerfile; end
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = []; end
if ~isfield(cfg, 'dataformat'),   	cfg.dataformat = []; end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');
cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile);

disp(hdr);

% open the figure and set the colormap
h = figure;
colormap(gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blocksize  = 1;
prevSample = 0;

while true
  
  % determine number of samples available in buffer
  hdr = read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);
  
  if hdr.nSamples > prevSample
    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample+1;
      endsample  = prevSample+blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end
    
    prevSample  = endsample;
    fprintf('\nprocessing scan %d\n', begsample);
    
    % read data from buffer (only the last scan)
    dat = read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample);
    
    % read events from buffer
    ev = read_event(cfg.eventfile, 'header', hdr);
    ev_sam = [ev.sample];
    ind = find(ev_sam==begsample-1);
    
    if ~isempty(ind)
      disp([ev(ind).type ev(ind).value])
    end
    
    % go to the correct figure
    figure(h);
    
    % this is a bit of a hack right now
    % information about number of slices etc. should go into the header
    % right now, pixeldata is transmitted as it's in the files
    W = sqrt(hdr.nChans);
    imagesc(reshape(dat,W,W)',[0 2048]);
    
    % force Matlab to update the figure
    drawnow
  else
    pause(0.01);
  end % if enough new samples
end % while true
