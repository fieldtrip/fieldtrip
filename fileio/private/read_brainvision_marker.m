function [rej] = read_brainvision_marker(fn);

% READ_BRAINVISION_MARKER reads rejection marks from a BrainVision file
%
% Use as
%   [rej] = read_brainvision_marker(filename)
%
% This function returns a Nx2 matrix with the begin and end latency
% of N rejection marks. The latency is in miliseconds.

% Copyright (C) 2004, Robert Oostenveld & Doug Davidson
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
% $Id: read_brainvision_marker.m 945 2010-04-21 17:41:20Z roboos $

rej = [];

[x y] = textread(fn, '%s%s', 1, 'delimiter', ',:' );

p   = mat2str(cell2mat(y));
pos = findstr('Hz', p);

samplingrate = str2num(p(1:pos-1));

s = 1./samplingrate;

[col1 col2] = textread(fn, '%*s%*s%n%n%*[^\n]',...
    'delimiter', ',',...
    'headerlines', 2 );

% Start and end points in msec
startp = (col1*s)*1000;
endp = ((col2*s)+(col1*s))*1000;

rej = [startp endp];
clear x p pos samplingrate;
