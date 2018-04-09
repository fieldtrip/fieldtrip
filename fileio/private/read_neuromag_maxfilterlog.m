function log = read_neuromag_maxfilterlog(filename)

% READ_NEUROMAG_MAXFILTERLOG reads the ascii logfile that is produced by MaxFilter
%
% Use as
%  log = read_neuromag_maxfilterlog(filename)
%
% See also READ_NEUROMAG_EVE, READ_NEUROMAG_HC

% Copyright (C) 2017, Robert Oostenveld
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
if fid<0
  ft_error('cannot open %s', filername);
end

% ************************************************************
% *** This is Elekta Neuromag MaxFilter(TM) version 2.2.15 ***
% *** Compilation date: Dec 11 2012 14:48:44               ***
% *** Please report problems to meg-support@elekta.com     ***
% ************************************************************

line = fgetl(fid);
line = fgetl(fid); assert(contains(line, 'MaxFilter'));
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);

numhpi = nan;

while ~feof(fid)
  % read the next line
  line = fgetl(fid);
  
  if startsWith(line, 'HPIFIT:') && contains(line, 'digitized')
    % HPIFIT: 4 coils digitized in order 1 2 3 4
    numhpi = sscanf(line, 'HPIFIT: %d');
    for i=1:numhpi
      log.hpi{i} = zeros(11, 0);
    end
    log.t = [];
    log.e = [];
    log.g = [];
    log.v = [];
    log.r = [];
    log.d = [];
  end
  
  if startsWith(line, 'Hpi fit OK')
    % Hpi fit OK, movements [mm/s] =  55.8 / 48.1 / 42.7 / 52.3 /
    %    79.5    0.0    0.0 /  79.2   -0.3    0.2 / g = 0.999 err =  0.5 mm / 0.998 -0.056 -0.011
    %     0.0   79.5    0.0 /   0.1   80.0   -0.2 / g = 1.000 err =  0.5 mm / 0.057  0.989  0.137
    %   -79.5    0.0    0.0 / -78.8   -1.0    0.2 / g = 0.994 err =  1.2 mm / 0.003 -0.138  0.990
    %     0.0  -79.5    0.0 /  -0.6  -78.6   -0.2 / g = 0.999 err =  1.0 mm /  -0.2    1.4    5.8
    % #t = 170.510, #e = 0.08 cm, #g = 0.998, #v = 0.75 cm/s, #r = 0.61 rad/s, #d = 0.26 cm
    for i=1:numhpi
      line = fgetl(fid);
      % sometimes the line starts with a '('
      if line(1)=='('
        value = sscanf(line, '( %f %f %f / %f %f %f / g = %f err = %f mm ) %f %f %f');
      else
        value = sscanf(line, '%f %f %f / %f %f %f / g = %f err = %f mm / %f %f %f');
      end
      log.hpi{i} = cat(2, log.hpi{i}, value);
    end
    line = fgetl(fid);
    value = sscanf(line, '#t = %f, #e = %f cm, #g = %f, #v = %f cm/s, #r = %f rad/s, #d = %f cm');
    log.t = cat(1, log.t, value(1));
    log.e = cat(1, log.e, value(2));
    log.g = cat(1, log.g, value(3));
    log.v = cat(1, log.v, value(4));
    log.r = cat(1, log.r, value(5));
    log.d = cat(1, log.d, value(6));
  end
  
end % while ~feof

fclose(fid);
