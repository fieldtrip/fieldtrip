function [dat] = read_labview_dtlg(filename, datatype)

% READ_LABVIEW_DTLG
%
% Use as
%   dat = read_labview_dtlg(filename, datatype)
% where datatype can be 'int32' or 'int16'
%
% The output of this function is a structure.

% Copyright (C) 2007, Robert Oostenveld
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


fid     = fopen(filename, 'r', 'ieee-be');

header  = fread(fid, 4, 'uint8=>char')';
if ~strcmp(header, 'DTLG')
  ft_error('unsupported file, header should start with DTLG');
end

version     = fread(fid, 4, 'char')'; % clear version
nd          = fread(fid, 1, 'int32');
p           = fread(fid, 1, 'int32');

% the following seems to work for files with version [7 0 128 0]
% but in files files with version [8 0 128 0] the length of the descriptor is not correct
ld          = fread(fid, 1, 'int16');
descriptor  = fread(fid, ld, 'uint8=>char')';

% ROBOOS: the descriptor should ideally be decoded, since it contains the variable
% name, type and size

% The first offset block always starts immediately after the data descriptor (at offset p, which should ideally be equal to 16+ld)
if nd<=128
  % The first offset block contains the offsets for all data sets.
  % In this case P points to the start of the offset block.
  fseek(fid, p, 'bof');
  offset = fread(fid, 128, 'uint32')';
else
  % The first offset block contains the offsets for the first 128 data sets.
  % The entries for the remaining data sets are stored in additional offset blocks.
  % The locations of those blocks are contained in a block table starting at P.
  offset = [];
  fseek(fid, p, 'bof');
  additional = fread(fid, 128, 'uint32');
  for i=1:sum(additional>0)
    fseek(fid, additional(i), 'bof');
    tmp    = fread(fid, 128, 'uint32')';
    offset = cat(2, offset, tmp);
  end
  clear additional i tmp
end

% ROBOOS: remove the zeros in the offset array for non-existing datasets
offset = offset(1:nd);

% ROBOOS: how to determine the data datatype?
switch datatype
  case 'uint32'
    datasize = 4;
  case 'int32'
    datasize = 4;
  case 'uint16'
    datasize = 2;
  case 'int16'
    datasize = 2;
  otherwise
    ft_error('unsupported datatype');
end

% If the data sets are n-dimensional arrays, the first n u32 longwords in each data
% set contain the array dimensions, imediately followed by the data values.

% ROBOOS: how to determine whether they are n-dimensional arrays?

% determine the number of dimensions by looking at the first array
% assume that all subsequent arrays have the same number of dimensions
if nd>1
  estimate = (offset(2)-offset(1)); % initial estimate for the number of datasize in the array
  fseek(fid, offset(1), 'bof');
  n = fread(fid, 1, 'int32');
  while mod(estimate-4*length(n), (datasize*prod(n)))~=0
    % determine the number and size of additional array dimensions
    n = cat(1, n, fread(fid, 1, 'int32'));
    if datasize*prod(n)>estimate
      ft_error('could not determine array size');
    end
  end
  ndim = length(n);
  clear estimate n
else
  estimate = filesize(fid)-offset;
  fseek(fid, offset(1), 'bof');
  n = fread(fid, 1, 'int32');
  while mod(estimate-4*length(n), (datasize*prod(n)))~=0
    % determine the number and size of additional array dimensions
    n = cat(1, n, fread(fid, 1, 'int32'));
    if datasize*prod(n)>estimate
      ft_error('could not determine array size');
    end
  end
  ndim = length(n);
  clear estimate n 
end

% read the dimensions and the data from each array
for i=1:nd
  fseek(fid, offset(i), 'bof');
  n = fread(fid, ndim, 'int32')';
  % Labview uses the C-convention for storing data, and MATLAB uses the Fortran convention
  n = fliplr(n);
  data{i} = fread(fid, n, datatype);
end
clear i n ndim

fclose(fid);

% put all local variables into a structure, this is a bit unusual programming style
% the output structure is messy, but contains all relevant information
tmp = whos;
dat = [];
for i=1:length(tmp)
  if isempty(strmatch(tmp(i).name, {'tmp', 'fid', 'ans', 'handles'}))
    dat = setfield(dat, tmp(i).name, eval(tmp(i).name));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local helper function to determine the size of the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function siz = filesize(fid)
fp = ftell(fid);
fseek(fid, 0, 'eof');
siz = ftell(fid);
fseek(fid, fp, 'bof');

