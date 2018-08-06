function [nvar, names, types] = nex_info(filename)
% nex_info(filename) -- read and display .nex file info
%
% [nvar, names, types] = nex_info(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% OUTPUT:
%   nvar - number of variables in the file
%   names - [nvar 64] array of variable names
%   types - [1 nvar] array of variable types
%           Interpretation of type values: 0-neuron, 1-event, 2-interval, 3-waveform, 
%                        4-population vector, 5-continuous variable, 6 - marker

% original from Plexon, download from http://www.plexoninc.com (8/4/02)
% modifications by Robert Oostenveld
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

if(nargin ~= 1)
   disp('1 input arguments are required')
   return
end

if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
    filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r', 'ieee-le');
if(fid == -1)
    disp('cannot open file');
   return
end

disp(strcat('file = ', filename));
magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
disp(strcat('version = ', num2str(version)));
disp(strcat('frequency = ', num2str(freq)));
disp(strcat('duration (sec) = ', num2str((tend - tbeg)/freq)));
disp(strcat('number of variables = ', num2str(nvar)));
names = zeros(1, 64);
for i=1:nvar
    types(i) = fread(fid, 1, 'int32');
    var_version = fread(fid, 1, 'int32');
    names(i, :) = fread(fid, [1 64], 'char');
    dummy = fread(fid, 128+8, 'char');
end
names = char(names);
fclose(fid);
