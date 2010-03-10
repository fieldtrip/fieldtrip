function fmri_viewer(cfg)

if ~isfield(cfg, 'headerfile'),    	cfg.headerfile = 'buffer://localhost:1972'; end
if ~isfield(cfg, 'datafile'),    	cfg.datafile   = cfg.headerfile; end
if ~isfield(cfg, 'eventfile'),    	cfg.eventfile  = cfg.headerfile; end
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = []; end
if ~isfield(cfg, 'dataformat'),   	cfg.dataformat = []; end

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile);

disp(hdr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prevSample = 0;

while true

  % determine number of samples available in buffer
  hdr = read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);

  if hdr.nSamples > prevSample
    % determine the samples to process
    begsample  = hdr.nSamples;
	endsample  = hdr.nSamples;

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
    	
    % this is a bit of a hack right now
    % information about number of slices etc. should go into the header
    % right now, pixeldata is transmitted as it's in the files
	W = sqrt(hdr.nChans);
    imagesc(reshape(dat,W,W)',[0 2048]);
	colormap(gray);

    % force Matlab to update the figure
	drawnow
  else
    pause(0.01);
  end % if enough new samples
end % while true
