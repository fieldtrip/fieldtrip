function fmri_viewer(cfg)

if ~isfield(cfg, 'headerfile'),    	cfg.headerfile = 'buffer://localhost:1972'; end
if ~isfield(cfg, 'datafile'),    	cfg.datafile = 'buffer://localhost:1972'; end
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
    fprintf('processing scan %d\n', begsample);

    % read data segment from buffer
    dat = read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample);

    % it only makes sense to read those events associated with the currently processed data
	%if strcmp(cfg.readevent, 'yes')
    %  evt = read_event(cfg.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);
    %end
	
	W = sqrt(hdr.nChans);
    imagesc(reshape(dat,W,W)',[0 2048]);
	colormap(gray);

    % force Matlab to update the figure
	drawnow

  end % if enough new samples
end % while true
