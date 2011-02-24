function [dat] = read_wdq_data(filename, hdr, begsample, endsample, chanindx)

% READ_WDQ_DATA reads data from wdq files
%
% Use as
%  [dat] = read_wdq_data(filename, hdr, begsample, endsample, chanindx)
%

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

fid = fopen(filename, 'r');

% set file pointer to where the data starts in the file
offset   = begsample * hdr.nchan * 2;
nsamples = endsample - begsample + 1;
fseek(fid, hdr.nbyteshdr + offset, 'bof');

datorig = fread(fid, [hdr.nchan nsamples], 'uint16');

% the lowest two bits may contain event info
lowbits = bitand(datorig, 1) + bitand(datorig, 2);
dat     = datorig - 2*bitand(datorig, 2^15);

% the higher 14 bits contain the waveforms
dat     = bitshift(dat, -2);

fclose(fid);

% scale the data
dat     = diag([hdr.chanhdr(1:hdr.nchan).scale]) * dat + ...
            [hdr.chanhdr(1:hdr.nchan).intercept]' * ones(1,size(dat,2));
          
dat     = dat(chanindx, :);