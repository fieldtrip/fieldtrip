function data = read_asa_msr(fn);

% READ_ASA_MSR reads EEG or MEG data from an ASA data file
% converting the units to uV or fT

% Copyright (C) 2002, Robert Oostenveld
% 
% $Log: read_asa_msr.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.5  2008/11/14 07:41:17  roboos
% only whitespace, no functional change
%
% Revision 1.4  2008/11/14 07:36:24  roboos
% use strcmpi instead of strcmp(lower())
%
% Revision 1.3  2008/11/12 17:02:03  roboos
% explicitely specify ieee-le in fopen()
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

Npnt      = read_asa(fn, 'NumberPositions=', '%d');
Ntime     = read_asa(fn, 'NumberTimesteps=', '%d');
UnitT     = read_asa(fn, 'UnitTime', '%s');
UnitM     = read_asa(fn, 'UnitMeas', '%s');
Timesteps = read_asa(fn, 'Timesteps', '%s');
lab       = read_asa(fn, 'Labels', '%s', Npnt);

val = read_asa(fn, 'Values', '%f');
if any(size(val)~=[Npnt,Ntime])
  msm_file = read_asa(fn, 'Values', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen(fullfile(path, msm_file), 'rb', 'ieee-le');
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
  error(sprintf('Unknown unit of time (%s)', UnitT));
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
  error(sprintf('Unknown unit of measurement (%s)', UnitM));
end

if length(size(lab))==2
  lab = tokenize(lab{1});
end

data.time  = time;
data.data  = val;
data.label = lab;


