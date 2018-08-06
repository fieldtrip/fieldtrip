function [el, lab] = read_elec(fn)

% READ_ELEC reads "la/mu" electrode parameters from a MBF electrode file
% which are used to position them on a triangulated surface
%
% [el, lab] = read_elec(filename)
%
% where el = [tri, la, mu]
% and lab contains the electrode labels (if present)
%
% See also READ_TRI, TRANSFER_ELEC

% Copyright (C) 1998, Robert Oostenveld
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

fid = fopen(fn, 'rt');
if fid~=-1

  % read the number of electrodes
  Nel = sscanf(fgetl(fid), '%d'); 
 
  % read the electrode triangle, lambda and mu
  for i=1:Nel
    str = fgetl(fid);
    el(i,:)  = sscanf(str, '%f %f %f')';
    indx = find(str=='!');
    if (indx)
      lab(i,:) = sprintf('%6s', str((indx+1):length(str)));
    end
  end
  fclose(fid);

else
  ft_error('unable to open file');
end

