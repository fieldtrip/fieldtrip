function data = read_asa_msr(fn)

% READ_ASA_MSR reads EEG or MEG data from an ASA data file
% converting the units to uV or fT

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

Npnt      = read_ini(fn, 'NumberPositions=', '%d');
Ntime     = read_ini(fn, 'NumberTimesteps=', '%d');
UnitT     = read_ini(fn, 'UnitTime', '%s');
UnitM     = read_ini(fn, 'UnitMeas', '%s');
Timesteps = read_ini(fn, 'Timesteps', '%s');
lab       = read_ini(fn, 'Labels', '%s', Npnt);

val = read_ini(fn, 'Values', '%f');
if any(size(val)~=[Npnt,Ntime])
  msm_file = read_ini(fn, 'Values', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen_or_error(fullfile(path, msm_file), 'rb', 'ieee-le');
  val = fread(fid, [Ntime, Npnt], 'float32')';
  fclose(fid);
end

tmp = sscanf(Timesteps, '%f(%f)%f');
time = linspace(tmp(1), tmp(3), Ntime);

if strcmpi(UnitT,'ms')
  time = 1*time;
elseif strcmpi(UnitT,'s')
  time = 1000*time;
elseif ~isempty(UnitT)
  ft_error(sprintf('Unknown unit of time (%s)', UnitT));
end

if strcmpi(UnitM,'uv')
  val = 1*val;
elseif strcmpi(UnitM,'?v')
  val = 1*val;
elseif strcmpi(UnitM,'mv')
  val = 1000*val;
elseif strcmpi(UnitM,'v')
  val = 1000000*val;
elseif strcmpi(UnitM,'ft')
  val = 1*val;
elseif strcmpi(UnitM,'pt')
  val = 1000*val;
elseif ~isempty(UnitM)
  ft_error(sprintf('Unknown unit of measurement (%s)', UnitM));
end

if length(size(lab))==2
  lab = tokenize(lab{1});
end

data.time  = time;
data.data  = val;
data.label = lab;


