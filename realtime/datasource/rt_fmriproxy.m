function rt_fmriproxy(cfg)
% RT_FMRIPROXY simulates an fMRI acquisition system by writing pixeldata in
% a 2s cycle. The pixeldata in this case is treated as a column vector with
% 122880 channels, but actually this represents 30 slices of 64x64 pixel
% each, one slice after another. The pixeldata is the same in each cycle 
% apart from added noise.
%
% Use as
%   rt_fmriproxy(cfg)
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2010, Stefan Klanke / Robert Oostenveld

% set the defaults
if isempty(cfg) | ~isfield(cfg.target, 'datafile'),    cfg.target.datafile = 'buffer://localhost:1972';  end
cfg.target.dataformat = [];    

% this should produce a variable called 'pixeldata'
load('demo_pixeldata');
% we don't catch errors here, since we rely on the contents of the file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a fieldtrip compatible header structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr = [];
hdr.Fs                 = 0.5;
hdr.nChans             = size(PixelDataInSlices,1);
hdr.nSamples           = 0;                                   
hdr.nSamplesPre        = 0;
hdr.nTrials            = 1;                           
hdr.label              = [];
hdr.siemensap          = AsciiProtocol;

endsample = 0;
stopwatch = tic;

% prepare an event structure that will go along with every sample/scan
ev = struct('type','timestamp','value','','offset',0,'duration',0,'sample',0);

while true
	% simulate acquisition of pixeldata
	t0 = toc(stopwatch);
	pause((endsample+1)/hdr.Fs - t0);

	pix = int16(double(PixelDataInSlices) + 200*rand(size(PixelDataInSlices)));
  
  unixtime  = etime(clock, [1970 1 1 0 0 0]);
  ev.value  = sprintf('%.6f',unixtime);
  ev.sample = endsample;     
  
	endsample = endsample + 1;
  
	fprintf('number of samples acquired = %i, sample time = %f, clock time = %f\n', endsample, endsample/hdr.Fs, toc(stopwatch));

	if endsample==1
    % flush the file, write the header and subsequently write the data segment
    ft_write_data(cfg.target.datafile, pix, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
  else
    % write the data segment
    ft_write_data(cfg.target.datafile, pix, 'append', true);
  end 
  
  ft_write_event(cfg.target.datafile, ev);
	% TODO: send an event along each with scan

	hdr.nSamples = endsample;

end % while again
