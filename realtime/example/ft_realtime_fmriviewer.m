function ft_realtime_fmriviewer(cfg)

% FT_REALTIME_FMRIVIEWER
%
% Use as
%   ft_realtime_fmriviewer(cfg)
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
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile);

disp(hdr);

if isfield(hdr,'nifti_1')
	ddim = double(hdr.nifti_1.dim);
    width  = ddim(1);
    height = ddim(2);
    numSlices = ddim(3);
elseif isfield(hdr,'siemensap') && isstruct(hdr.siemensap)
    width = siemensap.sKSpace.lBaseResolution;
    phaseFOV = siemensap.sSliceArray.asSlice{1}.dPhaseFOV;
    readoutFOV = siemensap.sSliceArray.asSlice{1}.dReadoutFOV;
    height = width * phaseFOV / readoutFOV;
    numSlices = siemensap.sSliceArray.lSize;
else
  warning('No protocol information found!')
  width = sqrt(hdr.nChans);
  height = width;
  numSlices = 1;
end
  

% open the figure and set the colormap
h = figure;
colormap(gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blocksize  = 1;
prevSample = 0;

% simple smoothing kernel
kern1d = [1 2 1];
kern = convn(kern1d'*kern1d,reshape(kern1d,1,1,3))
kern = kern/sum(kern(:));

while true
  
  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);
  
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
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample);
    
    % read events from buffer
    ev = ft_read_event(cfg.eventfile, 'header', hdr);
    if ~isempty(ev)
      ev_sam = [ev.sample];
      ind = find(ev_sam==begsample-1);
      if ~isempty(ind)
        disp([ev(ind).type ev(ind).value])
      end
    end
    
    % go to the correct figure
    figure(h);
    
    % actually, you should reshape to 3D (remove the *)
    S = reshape(dat, [width  height*numSlices]);
    
    V = reshape(dat, [width  height numSlices]);
    tic;
    Vn = convn(V,kern);
    toc
    Vn = Vn(2:end-1,2:end-1,2:end-1);
    size(Vn)
    S = reshape(Vn, [width height*numSlices]);

    
    imagesc(S(:,:)',[0 2048]);
    
    % force Matlab to update the figure
    drawnow
  else
    pause(0.01);
  end % if enough new samples
end % while true
