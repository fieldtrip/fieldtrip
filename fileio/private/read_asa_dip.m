function [pos, mom, ampl, time] = read_asa_dip(fn)

% READ_ASA_DIP reads the dipole position, moment and amplitude
% This importer is designed for fixed-dipole models and only supports 
% a limited number of the options that ASA has.
%
% Use as
%   [pos, mom, ampl, time] = read_asa_dip(filename)
%
% See also READ_ASA_VOL, READ_ASA_MRI

% Copyright (C) 2000, Robert Oostenveld
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

Npos = 0;
Ntim = 0;
pos  = [];
mom  = [];
ampl = [];
time = [];
fid1 = fopen_or_error(fn, 'r');

while ~Npos
  line = fgetl(fid1);
  if strncmp('NumberPositions=', line, length('NumberPositions='))
    Npos = sscanf(line, 'NumberPositions=%d');
  end
end
frewind(fid1);

while ~Ntim
  line = fgetl(fid1);
  if strncmp('NumberTimeSteps=', line, length('NumberTimeSteps='))
    Ntim = sscanf(line, 'NumberTimeSteps=%d');
  end
end
frewind(fid1);

while isempty(time)
  line = fgetl(fid1);
  if strncmp('TimeSteps', line, length('TimeSteps'))
    t = sscanf(line, 'TimeSteps %f(%f)%f');
        time = t(1):t(2):t(3);
  end
end
frewind(fid1);

if length(time)~=Ntim
  ft_error('incorrect timescale');
end

while isempty(pos)
  line = fgetl(fid1);
  if strncmp('Positions', line, length('Positions'))
    pos = zeros(3, Npos);
    pos(:) = fscanf(fid1, '%f');
    pos = pos';
  end
end
frewind(fid1);

while isempty(mom)
  line = fgetl(fid1);
  if strncmp('MomentsFixed', line, length('MomentsFixed'))
    mom = zeros(3, Npos);
    mom(:) = fscanf(fid1, '%f');
    mom = mom';
  end
end
frewind(fid1);

while isempty(ampl)
  line = fgetl(fid1);
  if strncmp('Magnitudes', line, length('Magnitudes'))
    ampl = zeros(Ntim, Npos);
    ampl(:) = fscanf(fid1, '%f');
    ampl = ampl';
  end
end
frewind(fid1);

fclose(fid1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old output format
%   pos = [x1 y1 z1 x2 y2 z2             ... ] and
%
%       mom = [Qx11 Qy11 Qz11 Qx21 Qy21 Qz21 ...
%              Qx12 Qy12 Qz12 Qx22 Qy22 Qz22 ...
%              ...
%              Qx1t Qy1t Qz1t Qx2t Qy2t Qz2t ... ]
%
%   time = [T1 T2 T3 T4 .. Tt]
%
% was obtained by re-format the data according to
%
% pos = pos';
% pos = pos(:)';
% 
% mom1 = [];
% for j=1:Npos
%   mom1 = [mom1 ampl(j,:)' * mom(j,:)];
% end
% mom = mom1;

