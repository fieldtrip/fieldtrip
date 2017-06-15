function ft_realtime_fmriproxy(cfg)

% FT_REALTIME_FMRIPROXY simulates an fMRI acquisition system by writing volumes in a
% cycle of about 2 seconds. The voxel data is written as a column vector with X*Y*Z
% channels, where X and Y are the readout and phase-encoding resolution, and Z is the
% number of slices. The voxel data consists of a ellipsoid (a bit like a head) with
% added lateralized activation (20 scan cycle) and noise.
%
% This function also writes out events of type='scan' and value='pulse' when the
% simulated scanner initiates a scan, and value='ready' when a hypothetical
% processing pipeline is finished with that scan, just after writing out the volume
% data itself. There is an artificial delay of 1.3*TR between the two events.
%
% Use as
%   ft_realtime_fmriproxy(cfg)
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%
% You can also select a resolution of the simulated volumes (default = [64,64,20]) like
%   cfg.voxels = [64 64 32]
% and the repetition time (TR, default = 0.08*number of slices) in seconds using
%   cfg.TR = 2.0
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2010, Stefan Klanke / Robert Oostenveld
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
if isempty(cfg) || ~isfield(cfg, 'target') || ~isfield(cfg.target, 'datafile')
  cfg.target.datafile = 'buffer://localhost:1972';  
end
cfg.target.dataformat = [];    

if isfield(cfg, 'voxels')
   voxels = cfg.voxels(1:3);
else
   voxels = [64 64 20];
end

if isfield(cfg,'TR')
   TR = cfg.TR;
   Tslice = TR / voxels(3);
else
   Tslice = 0.08;
   TR = Tslice * voxels(3);
end

nifti = [];
nifti.dim = voxels;
nifti.pixdim = [3.5 3.5 3.5]; % size of voxel in mm
nifti.slice_duration = Tslice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a FieldTrip compatible header structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdr = [];
hdr.nChans = prod(voxels);
hdr.nSamples = 0;
hdr.Fs = 1/TR;
hdr.nSamplesPre        = 0;
hdr.nTrials            = 1;                           
hdr.nifti_1 = encode_nifti1(nifti);

% prepare two event structure that will go along with every sample/scan
evp = [];
evp.type = 'scan';
evp.value = 'pulse';
evp.offset = 0;
evp.duration = 0;
evp.sample = 0;

evr = evp;
evr.value = 'ready';

% assume acquisition + processing takes 30% of a TR cycle
TP = 0.3 * TR;

% make some fake MRI data
% Sphere = volume containing a spherical "head"
% Left/Right = volume of same size, adds a bit of fake activation
x = linspace(-1,1,voxels(1));
y = linspace(-1,1,voxels(2));
z = linspace(-1,1,voxels(3));
[X,Y,Z] = meshgrid(x,y,z);

R = sqrt(X.^2 + Y.^2 + Z.^2);
Sphere = single(1200 * (1./(1 + exp(40*(R-0.8)))));
R = sqrt((X-0.4).^2 + Y.^2 + Z.^2);
Left  = single(60 * (1./(1 + exp(40*(R-0.4)))));
R = sqrt((X+0.4).^2 + Y.^2 + Z.^2);
Right = single(60 * (1./(1 + exp(40*(R-0.4)))));

stopwatch = tic;

% write header to FieldTrip buffer
ft_write_data(cfg.target.datafile, single([]), 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);

numPulse = 0;
numScan = 0;

while true
	% simulate acquisition of pixeldata
	t0 = toc(stopwatch);
	
	% write 'scan pulse' event
	numPulse = numPulse + 1;
	evp.sample = numPulse;
	fprintf('Sending pulse: %3i, clock time = %f\n', numPulse, t0);	
	ft_write_event(cfg.target.datafile, evp);
	
	if numPulse > 1
		% fake acquisition + processing
		lateral = 0.5 + 0.5*sin(2*pi*numScan/20);
	
		Noise = single(30*randn(voxels));
		scan = single(Sphere + lateral*Left + (1-lateral)*Right + Noise);
	
		tp = toc(stopwatch);
		pause(TP - (tp-t0));
		
		ft_write_data(cfg.target.datafile, scan(:), 'append', true);
		% write 'scan ready' event
		numScan = numScan + 1;
		evr.sample = numScan;
		fprintf('Scans written: %3i, clock time = %f\n', numScan, toc(stopwatch));	
		ft_write_event(cfg.target.datafile, evr);
	end
	
	tr = toc(stopwatch);
	pause(numPulse*TR - tr);
end
