function rt_fmriproxy(cfg)
% RT_FMRIPROXY simulates an fMRI acquisition system by writing pixeldata in
% a 2s cycle. The pixeldata in this case is treated as a column vector with
% 147456 channels, but actually this represent a 6x6 tiled image of 64x64 pixel
% slices. The pixeldata is the same in each cycle apart from added noise.
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

try
  % this should produce a variable called 'pixeldata'
  load('demo_pixeldata');
catch
end

% file not found? just create noise in the right scale
if ~exist('pixeldata', 'var')
	pixeldata = 1500*rand(147456,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a fieldtrip compatible header structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr = [];
hdr.Fs                 = 0.5;
hdr.nChans             = size(pixeldata,1);
hdr.nSamples           = 0;                                   
hdr.nSamplesPre        = 0;
hdr.nTrials            = 1;                           
hdr.label              = [];

endsample = 0;
stopwatch = tic;

while true
	% simulate acquisition of pixeldata
	t0 = toc(stopwatch);
	pause((endsample+1)/hdr.Fs - t0);

	pix = int16(pixeldata + 200*rand(size(pixeldata)));
	endsample = endsample + 1;
	
	fprintf('number of samples acquired = %i, sample time = %f, clock time = %f\n', endsample, endsample/hdr.Fs, toc(stopwatch));

	if endsample==1
      % flush the file, write the header and subsequently write the data segment
      write_data(cfg.target.datafile, pix, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
    else
      % write the data segment
	  write_data(cfg.target.datafile, pix, 'append', true);
    end 
	% TODO: send an event along each with scan

	hdr.nSamples = endsample;

end % while again
