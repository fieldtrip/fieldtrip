function ft_realtime_pooraudioproxy(cfg)

% FT_REALTIME_POORAUDIOPROXY reads continuous data from the sound card using the
% standard Matlab API and writes it to a FieldTrip buffer. This proxy has poor timing
% and will produce dropped audio frames between blocks. Also the Matlab documentation
% warns about using this API for long recordings because this will fill up memory and
% degrade performance.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_pooraudioproxy(cfg)
%
% The audio-specific configuration structure can contain
%   cfg.channel     = number of channels (1 or 2, default=2)
%   cfg.blocksize   = size of recorded audio blocks in seconds (default=1)
%   cfg.fsample     = audio sampling frequency in Hz (default = 44100)
%   cfg.nbits       = recording depth in bits (default = 16)
%
% Note that currently, the sound will be buffered in double precision irrespective of the sampling bit depth.
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% Finally, there is an option for showing debug output
%   cfg.debug       = show sample time and clock time (default = 'yes')
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2010, Stefan Klanke & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the defaults
if ~isfield(cfg, 'target'),             cfg.target = [];                                  end
if ~isfield(cfg, 'blocksize'),          cfg.blocksize = 1;                                end % in seconds
if ~isfield(cfg, 'channel'),            cfg.channel = 2;                                  end % default is stereo
if ~isfield(cfg, 'fsample'),            cfg.fsample = 44100;                              end % in Hz
if ~isfield(cfg, 'nbits'),              cfg.nbits = 16;                                   end % default 16 bit
if ~isfield(cfg, 'debug'),              cfg.debug = 'yes';                                end
if ~isfield(cfg.target, 'datafile'),    cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg.target, 'dataformat'),  cfg.target.dataformat = [];                       end % default is to use autodetection of the output format

% it seems that audiorecorder can only record integer multiples of 0.5 seconds
% this was observed on OSX with Matlab 2008a, but has not been confirmed on
% other platforms
newblocksize = round(cfg.blocksize*2)/2; % round to nearest half or full second
if newblocksize<0.5
  newblocksize = 0.5;
end
if newblocksize~=cfg.blocksize
  warning('sestting cfg.blocksize to %f', newblocksize);
  cfg.blocksize = newblocksize;
end

hdr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a FieldTrip compatible header structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr.Fs                 = cfg.fsample;								  
hdr.nChans             = cfg.channel;					                  
hdr.nSamples           = 0;                                   
hdr.nSamplesPre        = 0;
hdr.nTrials            = 1;                           
if cfg.channel == 2
	hdr.label              = {'Audio Left','Audio Right'};
elseif cfg.channel == 1
	hdr.label              = {'Audio Left'};
else
	error 'Invalid channel number (only 1 or 2 are allowed)';
end
hdr.FirstTimeStamp     = nan;
hdr.TimeStampPerSample = nan;

% create an audio recorder object
REC = audiorecorder(cfg.fsample, cfg.nbits, cfg.channel);

count = 0;
endsample = 0;
stopwatch  = tic;

while true
	% tell the audio recorder to read (and wait for) cfg.blocksize samples, and then get the data into the matrix x
	recordblocking(REC, cfg.blocksize);	
	x = getaudiodata(REC);
	endsample = endsample + size(x,1);
	
	% since this is a "poor" audio proxy, the clock time will move faster than the sample time here
	% in general, those two should be the same (or have a small constant offset)
	if strcmp(cfg.debug, 'yes')
		fprintf('number of samples acquired = %i, sample time = %f, clock time = %f\n', endsample, endsample/hdr.Fs, toc(stopwatch));
	end

	% increase the internal counter of blocks, so we know whether to append the new data or not
	count = count+1;
	
	if count==1
    % flush the file, write the header and subsequently write the data segment
    ft_write_data(cfg.target.datafile, x', 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
  else
    % write the data segment
    ft_write_data(cfg.target.datafile, x', 'append', true);
  end % if count==1

	hdr.nSamples = endsample;
	
end % while again
