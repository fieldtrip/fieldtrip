function [t] = neuralynx_timestamp(filename, num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading a single timestamp of a single channel Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

headersize = 16384;
switch ft_filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nst'
    recordsize = 304;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
  case 'neuralynx_ntt'
    recordsize = 304;   % in bytes
end

fid = fopen(filename, 'rb', 'ieee-le');

if (ispc)
  % this is to fix a bug in the windwos version which does not want to do uint64=>uint64
  % however this code will fail if the MSB is set (only likely in very long recordings)
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>integer*8');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>integer*8');
  end
  t = uint64(t);
else
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>uint64');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>uint64');
  end
end

fclose(fid);
