function [avg] = read_ns_avg(filename)

% READ_NS_AVG read a NeuroScan 3.x or 4.x AVG File
%
% [avg] = read_ns_avg(filename)
%
%  The output data structure avg has the fields:
%   avg.data        - ERP signal in uV (Nchan x Npnt)
%   avg.nsweeps     - number of accepted trials/sweeps in avg
%   avg.variance    - variance of the signal (Nchan x Npnt)
%   avg.label       - electrode labels
%   avg.nchan       - number of channels
%   avg.npnt        - number of samplepoints in ERP waveform
%   avg.rate        - sample rate (Hz)
%   avg.time        - time for each sample OR
%   avg.frequency   - frequency for each sample
%   hdr.domain      - flag indicating time (0) or frequency (1) domain
%   avg.xmin        - prestimulus epoch start (e.g., -100 msec)
%   avg.xmax        - poststimulus epoch end (e.g., 900 msec)

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_ns_avg.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2004/06/28 07:36:53  roberto
% changed from DOS to UNIX linefeeds, I am not completely sure whether I made other changes as well
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

% read the neuroscan header 
avg = read_ns_hdr(filename);

% create a time or frequency axis
if avg.domain==1
  avg.frequency = linspace(avg.xmin, avg.xmax, avg.npnt) / 1000;	% in Hz instead of mili-Hz
else
  avg.time = linspace(avg.xmin, avg.xmax, avg.npnt);			% in ms
end

% open the file and seek towards the place where the raw data is
fid = fopen(filename,'r','ieee-le');
if fid<0
  error(['cannot open ', filename]);
else
  fseek(fid, 900, 'cof');		% skip general header
  fseek(fid, 75*avg.nchan, 'cof');	% skip channel headers
end;

% read raw signal data and convert to uV
avg.data = zeros(avg.nchan, avg.npnt);
for elec = 1:avg.nchan
    fseek(fid, 5, 'cof'); % skip sweeps header
    raw = fread(fid, avg.npnt, 'float32');
    avg.data(elec,:) = (raw' - avg.baseline(elec)) * avg.calib(elec) / avg.nsweeps;
end;

% read signal variance if present
if avg.variance
    variance = zeros(avg.npnt, avg.nchan);
    for elec = 1:avg.nchan,
        variance(:, elec) = fread(fid, avg.npnt, 'float32');
    end;
    avg.variance = variance';
else
    avg.variance = [];
end;

fclose(fid);

