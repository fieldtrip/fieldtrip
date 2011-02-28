function [dat] = read_wdq_data(filename, hdr, begsample, endsample, chanindx)

% READ_WDQ_DATA reads data from wdq files
%
% Use as
%  [dat] = read_wdq_data(filename, hdr, begsample, endsample, chanindx)
% or
%  [dat] = read_wdq_data(filename, hdr, 'lowbits')

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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
% $Id: read_wdq_data.m $

% information about how to interpret the file are taken from the document
% 'CODAS data storage format'

if nargin==3 && ischar(begsample) && strcmp(begsample, 'lowbits')
  % read in the whole stretch of data and return the lowest 2 bits
  begsample = 1;
  endsample = 0.5*hdr.nbytesdat/hdr.nchan;
  getlowbits = 1;
elseif nargin==5
  getlowbits = 0;
end

fid = fopen(filename, 'r');

% set file pointer to where the data starts in the file
offset   = begsample * hdr.nchan * 2;
nsamples = endsample - begsample + 1;
fseek(fid, hdr.nbyteshdr + offset, 'bof');

datorig = fread(fid, [hdr.nchan nsamples], 'uint16');
fclose(fid);

% the lowest two bits may contain event info
% the higher 14 bits contain the waveforms
if getlowbits
  dat = bitand(datorig, 1) + bitand(datorig, 2);
else
  dat = datorig - 2*bitand(datorig, 2^15);
  dat = bitshift(dat, -2);

  % scale the data
  dat     = diag([hdr.chanhdr(1:hdr.nchan).scale]) * dat + ...
              [hdr.chanhdr(1:hdr.nchan).intercept]' * ones(1,size(dat,2));          
  dat     = dat(chanindx, :);
end