function [condNumbers,condLabels] = read_ctf_cls(fname)

% READ_CTF_CLS reads the classification file from a CTF dataset

% Copyright (C) 2003, Ole Jensen
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

condNumbers = [];
condLabels = {};

try
  fid = fopen_or_error(fname,'r');
catch err
  condNumbers = [];
  condLabels = [];
  return
end

nCondition = 0;
readBad = 0;
readList = 0;
S2 = '*';
S1 = '*';
while ~isempty(S1)
  S3 = S2;
  S2 = S1;
  S1 =fscanf(fid,'%s',1);

  if readList
    if ~isempty(S1) && ~isempty(str2num(S1(2:end)))
      k = k + 1;
      condTmp = [condTmp 1+str2num(S1(2:end))];
    else
      readList = 0;
    end
    condNumbers{nCondition} = condTmp;
  end

  if strcmp(S2,'NAME:')
    % New condition found!
    % fprintf('%s\n',S1);
    nCondition = nCondition+1;
    condLabels(nCondition) = {S1} ;
  end

  if strcmp(S1,'NUMBER') && strcmp(S2,'TRIAL')
    if ~isempty(S1)
      readList = 1;
      k = 0;
      condTmp = [];
    else
      readList = 0;
    end
  end

end % ~isempty(S1)

