function [t] = neuralynx_numrecords(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for determining the number of records in a single channel Neualynx file
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
% $Id: neuralynx_numrecords.m 945 2010-04-21 17:41:20Z roboos $

headersize = 16384;
switch filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
  otherwise
    error('unsupported filetype');
end

fid = fopen(filename, 'rb', 'ieee-le');
fseek(fid, 0, 'eof');
t = floor((ftell(fid) - headersize)/recordsize);
fclose(fid);
