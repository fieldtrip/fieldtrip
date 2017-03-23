function [magn] = read_ctf_sens(filename)

% READ_CTF_SENS reads MEG sensor information from CTF configuration file
%
% magn = read_ctf_sens(filename)
%
% where the returned structure meg has the fields
%   magn.pnt    position first coil
%   magn.ori    orientation first coil
%   magn.pnt2   position second coil
%   magn.ori2   orientation second coil

% Copyright (C) 2002, Robert Oostenveld
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

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

% skip the first line
line = fgetl(fid);

% read the channel parameters
chan = 0;
while (1)
  line = fgetl(fid);
  if ~isempty(line) & line==-1
    % reached end of file
    break
  end

  element = tokenize(line);
  if ~isempty(element)
    chan = chan+1;

    % channel label or name
    name{chan} = strtok(element{1}, '-');

    % type of channel (EEG/MEG etc.), convert to uppercase
    type{chan} = upper(element{7});

    % number of coils (magnetometer=1/gradiometer=2)
    ncoils(chan) = str2num(element{2});

    % polarity of channel
    polarity(chan) = str2num(element{26});

    % position of first coil
    coil1.pnt(chan,:) = [str2num(element{17}) str2num(element{18}) str2num(element{19})];
    % direction of first coil
    coil1.dir(chan,:) = [str2num(element{20}) str2num(element{21})];

    % length of basline
    baseline.length(chan) = str2num(element{16});
    % direction of baseline
    baseline.dir(chan,:) = [str2num(element{22}) str2num(element{23})];

    % direction of second coil
    coil2.dir(chan,:) = [str2num(element{24}) str2num(element{25})];
  end
end

% find the sensor gradiometer channels
indx = find(strcmp(type, 'GRAD1-SENS'));
nchan = length(indx);
label = name(indx);
pnt1 = coil1.pnt(indx,:) * 10;  % convert from cm to mm

% shift offset to origin to make triangulation
offset = [mean(pnt1(:,1)) mean(pnt1(:,2)) min(pnt1(:,3))];
prj = elproj(pnt1 - repmat(offset, nchan, 1));
tri = delaunay(prj(:,1), prj(:,2));

% compute the directional vectors of each gradiometer: bottom coil
a = coil1.dir(indx,1);
b = coil1.dir(indx,2);
[r1, r2, r3] = sph2cart(a*2*pi/360, pi/2 - b*2*pi/360, ones(nchan,1));
r1 = -r1 .* polarity(indx)';
r2 = -r2 .* polarity(indx)';
r3 = -r3 .* polarity(indx)';
ori1 = [r1 r2 r3];

% compute the directional vectors of each gradiometer: top coil
a = coil2.dir(indx,1);
b = coil2.dir(indx,2);
[r1, r2, r3] = sph2cart(a*2*pi/360, pi/2 - b*2*pi/360, ones(nchan,1));
r1 = -r1 .* polarity(indx)';
r2 = -r2 .* polarity(indx)';
r3 = -r3 .* polarity(indx)';
ori2 = [r1 r2 r3];

% compute the directional vectors of the baseline and the position of the 2nd coil
a = baseline.dir(indx,1);
b = baseline.dir(indx,2);
[r1, r2, r3] = sph2cart(a*2*pi/360, pi/2 - b*2*pi/360, ones(nchan,1));
oriB = [r1 r2 r3];
pnt2 = pnt1 + oriB .* repmat(baseline.length(indx)', 1, 3);

magn.pnt    = pnt1;
magn.ori    = ori1;
magn.pnt2   = pnt2;
magn.ori2   = ori2;
magn.tri    = tri;
magn.label  = label;
