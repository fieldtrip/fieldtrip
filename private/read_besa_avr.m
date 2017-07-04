function [avr] = read_besa_avr(filename)

% READ_BESA_AVR reads average EEG data in BESA format
%
% Use as
%   [avr] = read_besa_avr(filename)
%
% This will return a structure with the header information in
%   avr.npnt
%   avr.tsb
%   avr.di
%   avr.sb
%   avr.sc
%   avr.Nchan   (optional)
%   avr.label   (optional)
% and the ERP data is contained in the Nchan X Nsamples matrix
%   avr.data

% Copyright (C) 2003-2006, Robert Oostenveld
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

fid = fopen(filename, 'rt');

% the first line contains header information
headstr = fgetl(fid);
ok = 0;

if ~ok
  try
    buf = sscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f Nchan= %f SegmentName= %s\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    avr.Nchan = buf(6);
    avr.SegmentName = buf(7);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  try
    buf = fscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f Nchan= %f\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    avr.Nchan = buf(6);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  try
    buf = sscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  ft_error('Could not interpret the header information.');
end

% rewind to the beginning of the file, skip the header line
fseek(fid, 0, 'bof');
fgetl(fid);

% the second line may contain channel names
chanstr = fgetl(fid);
chanstr = deblank(fliplr(deblank(fliplr(chanstr))));
if (chanstr(1)>='A' && chanstr(1)<='Z') || (chanstr(1)>='a' && chanstr(1)<='z')
  haschan   = 1;
  avr.label = str2cell(strrep(deblank(chanstr), '''', ''))';
else
  [root, name] = fileparts(filename);
  haschan = 0;
  elpfile = fullfile(root, [name '.elp']);
  elafile = fullfile(root, [name '.ela']);
  if exist(elpfile, 'file')
    % read the channel names from the accompanying ELP file
    lbl = importdata(elpfile);
    avr.label = strrep(lbl.textdata(:,2) ,'''', '');
  elseif exist(elafile, 'file')
    % read the channel names from the accompanying ELA file
    lbl = importdata(elafile);
    lbl = strrep(lbl ,'MEG ', ''); % remove the channel type
    lbl = strrep(lbl ,'EEG ', ''); % remove the channel type
    avr.label = lbl;
  else
    ft_warning('Could not create channels labels.');
  end
end

% seek to the beginning of the data
fseek(fid, 0, 'bof');
fgetl(fid);   % skip the header line
if haschan
  fgetl(fid); % skip the channel name line
end

buf = fscanf(fid, '%f');
nchan = length(buf)/avr.npnt;
avr.data = reshape(buf, avr.npnt, nchan)';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to cut a string into pieces at the spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = str2cell(s)
c = {};
[t, r] = strtok(s, ' ');
while ~isempty(t)
  c{end+1} = t;
  [t, r] = strtok(r, ' ');
end
