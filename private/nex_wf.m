function [adfreq, n, ts, nf, w] = nex_wf(filename, varname)
% nex_wf(filename, varname): Read waveform variable from a .nex file
%
% [adfreq, n, ts, nf, w] = nex_wf(filename, varname)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
%
%           
% OUTPUT:
%   n - number of waveforms
%   ts - array of waveform timestamps (in seconds)
%   nf - number of data points in each waveform
%   w - matrix of waveform a/d values [n nf] (in millivolts)

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

n = 0;
adfreq = 0;
ts = 0;
nf = 0;
w = 0;

if(nargin ~= 2)
   disp('2 input arguments are required')
   return
end

if(ischar(filename) == 0)
   disp('input arguments should be character arrays')
   return
end

if(ischar(varname) == 0)
   disp('input arguments should be character arrays')
   return
end

if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
    filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r', 'ieee-le');
if(fid == 0)
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
name = zeros(1, 64);
found = 0;
for i=1:nvar
    type = fread(fid, 1, 'int32');
    var_version = fread(fid, 1, 'int32');
    name = fread(fid, [1 64], 'char');
    offset = fread(fid, 1, 'int32');
    nf = fread(fid, 1, 'int32');
    dummy = fread(fid, 32, 'char');
    adfreq = fread(fid, 1, 'double');
    adtomv = fread(fid, 1, 'double');
    n = fread(fid, 1, 'int32');
    name = char(name);
    name = deblank(name);
    k = strcmp(name, deblank(varname));
    if(k == 1)
        if type ~= 3
            disp(sprintf('%s is not a waveform variable', deblank(varname)));
            return;
        end
        found = 1;
        fseek(fid, offset, 'bof');
        ts = fread(fid, [1 nf], 'int32');
        w = fread(fid, [n nf], 'int16');
        break
    end
    dummy = fread(fid, 76, 'char');
end

fclose(fid);

if found == 0
    disp('did not find variable in the file');
else
    ts = ts/freq;
    w = w*adtomv;
    disp(strcat('number of waveforms = ', num2str(nf)));
end
