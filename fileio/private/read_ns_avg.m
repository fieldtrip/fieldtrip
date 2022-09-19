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

% Copyright (C) 2002, Robert Oostenveld
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

% read the neuroscan header 
avg = read_ns_hdr(filename);

% create a time or frequency axis
if avg.domain==1
  avg.frequency = linspace(avg.xmin, avg.xmax, avg.npnt) / 1000;    % in Hz instead of mili-Hz
else
  avg.time = linspace(avg.xmin, avg.xmax, avg.npnt);            % in ms
end

% open the file and seek towards the place where the raw data is
fid = fopen_or_error(filename,'r','ieee-le');
fseek(fid, 900, 'cof');       % skip general header
fseek(fid, 75*avg.nchan, 'cof');  % skip channel headers

% read raw signal data and convert to uV
avg.data = zeros(avg.nchan, avg.npnt);
for elec = 1:avg.nchan
    fseek(fid, 5, 'cof'); % skip sweeps header
    raw = fread(fid, avg.npnt, 'float32');
    avg.data(elec,:) = (raw' - avg.baseline(elec)) * avg.calib(elec) / avg.nsweeps;
end

% read signal variance if present
if avg.variance
    variance = zeros(avg.npnt, avg.nchan);
    for elec = 1:avg.nchan,
        variance(:, elec) = fread(fid, avg.npnt, 'float32');
    end
    avg.variance = variance';
else
    avg.variance = [];
end

fclose(fid);

