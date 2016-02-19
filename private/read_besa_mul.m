function [dat] = read_besa_mul(filename)

% READ_BESA_MUL reads data from a BESA multiplexed (*.mul) file
%
% Use as
%   dat = read_besa_mul(filename);

% Copyright (C) 2005, Robert Oostenveld
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

dat = [];
fid = fopen(filename, 'rt');

% According to the BESA documentation, the ASCII Multiplexed Format is as
% follows:
%
% The first of two header lines contains similar information to that of the
% BESA ASCII file:
%
% TimePoints= 200 Channels= 27 BeginSweep[ms]= -500.00
% SamplingInterval[ms]= 5.000 Bins/uV= 1.000 SegmentName=Condition1
%
% Note that the item 'SegmentName' is missing if no segment comment is
% specified when writing a segment to file.
%
% If an epoch of a continuous EEG is exported in ASCII multiplexed format,
% the first header line contains the additional item 'Time', which
% indicates the daytime of the first sample in the exported segment:
%
% TimePoints= 200 Channels= 27 BeginSweep[ms]= 0.00 SamplingInterval[ms]=
% 5.000 Bins/uV= 1.000 Time=22:02:53 SegmentName=Segment1

hdr1 = fgetl(fid);
% split the first header line into separate elements
%!!! 2009/09/25
%tmp = tokenize(hdr1, ' ');
tmp = tokenize(hdr1, ' ',1);
for i=1:length(tmp)
  % extract the information from each element
  dum = tokenize(tmp{i}, '=');
  var = dum{1}; % this is the name of the header variable
  val = dum{2}; % this is the value of the header variable
  var(var=='/') = '_';  % replace characters that are not valid in a structure element name
  var(var=='[') = '_';  % replace characters that are not valid in a structure element name
  var(var==']') = '_';  % replace characters that are not valid in a structure element name
  var(var=='.') = '_';  % replace characters that are not valid in a structure element name
  num = str2num(val);
  if ~isempty(num)
    dat = setfield(dat, var, num);
  else
    dat = setfield(dat, var, val);
  end
end

% The second line of the header contains labels for each channel, which may
% be either the original channel names, or the names of the channels of the
% current montage, e.g.
%
% O1 Oz P3 T5 T3 C3 F7 F3 Fp1 Fz Cz Pz Fp2 F4 F8 C4 T4 T6 P4 Fpz O2 M2 M1
% F10 F9 T10 T9
%
% Each subsequent line contains values for all 'Channels' at one time
% point, in floating point or scientific format. Values are given for the
% current or the original montage, selected as described above.
%
% Labels for source montages have the following form: 'TAr-L'.
%
% The first two letters indicate the head region:
%
% The small letter indicates in part the orientation: r=radial,
% t=tangential, and in part the relative location of the basal temporal
% source: l=lateral, m=mesial.
%
% The final letter after the hyphen indicates L=left, M=middle, R=right.

hdr2 = fgetl(fid);
% split the second header line into channel/source labels
%!!! 2009/09/25
%dat.label = tokenize(hdr2, ' ');
dat.label = tokenize(hdr2, ' ',1);

% read the actual data
dat.data = fscanf(fid, '%g');
dat.data = reshape(dat.data, [dat.Channels dat.TimePoints]);

fclose(fid);
return
