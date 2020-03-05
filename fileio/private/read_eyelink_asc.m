function asc = read_eyelink_asc(filename)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file
%
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010-2015, Robert Oostenveld
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

asc.header  = {};
asc.msg     = {};
asc.input   = [];
asc.sfix    = {};
asc.efix    = {};
asc.ssacc   = {};
asc.esacc   = {};
asc.dat     = [];
current   = 0;

% read the whole file at once
fid = fopen_or_error(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(newline));        % split on newline

for i=1:numel(aline)
  tline = aline{i};

  if numel(tline) && any(tline(1)=='0':'9')
  % if regexp(tline, '^[0-9]')
    tline   = strrep(tline, ' . ', ' NaN '); % replace missing values
    tmp     = sscanf(tline, '%f');
    nchan   = numel(tmp);
    current = current + 1;

    if size(asc.dat, 1)<nchan
      % increase the allocated number of channels
      asc.dat(nchan,:) = 0;
    end

    if size(asc.dat, 2)<current
      % increase the allocated number of samples
      asc.dat(:,end+10000) = 0;
    end

    % add the current sample to the data matrix
    asc.dat(1:nchan, current) = tmp;


  elseif regexp(tline, '^INPUT')
    [val, num] = sscanf(tline, 'INPUT %d %d');
    this.timestamp = val(1);
    this.value     = val(2);
    if isempty(asc.input)
      asc.input = this;
    else
      asc.input = cat(1, asc.input, this);
    end


  elseif regexp(tline, '\*\*.*')
    asc.header = cat(1, asc.header, {tline});


  elseif regexp(tline, '^MSG')
    asc.msg = cat(1, asc.msg, {tline});


  elseif regexp(tline, '^SFIX')
    asc.sfix = cat(1, asc.sfix, {tline});


  elseif regexp(tline, '^EFIX')
    asc.efix = cat(1, asc.efix, {tline});


  elseif regexp(tline, '^SSACC')
    asc.ssacc = cat(1, asc.ssacc, {tline});


  elseif regexp(tline, '^ESACC')
    asc.esacc = cat(1, asc.esacc, {tline});

  else
    % all other lines are not parsed
  end

end

% remove the samples that were not filled with real data
asc.dat = asc.dat(:,1:current);

