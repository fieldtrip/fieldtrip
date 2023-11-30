function [output] = volumeselectlargest(input)

% VOLUMESELECTLARGEST is a helper function for segmentations
%
% See also VOLUMEFILLHOLES, VOLUMETHRESHOLD, VOLUMESMOOTH, VOLUMEPAD

% Copyrights (C) 2023, Robert Oostenveld
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

ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

inflate = volumepad(input, 1);                    % pad the volume
[lab, num] = spm_bwlabel(double(inflate), 18);    % note that 18 is consistent with imfill, 26 is not

if num>1
  % determine the size of each cluster
  s = zeros(1, num);
  for i=1:num
    s(i) = sum(lab(:)==i);
  end

  % select the largest
  [maxsize, maxindex] = max(s);

  inflate = (lab==maxindex);
  output  = inflate(2:end-1, 2:end-1, 2:end-1);   % trim the edges
  
else
  output = input;
end