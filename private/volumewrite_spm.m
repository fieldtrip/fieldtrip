function [Va] = volumewrite_spm(filename, data, transform)

% VOLUMEWRITE_SPM write a anatomical or functional volume to img/hdr file
% using the SPM toolbox
%
% Use as
%   [Va] = volumewrite_spm(filename, data, transform)

% Copyright (C) 2006, Robert Oostenveld
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
% $Id$

% check whether the required SPM2 toolbox is available
hastoolbox('SPM2', 1);

if 0
  % not all datatypes are supported by SPM, this can be checked with
  spm_type('double') %           -- Double precision floating point number array
  spm_type('single') %           -- Single precision floating point number array
  spm_type('logical') %          -- Logical array
  spm_type('char') %             -- Character array
  spm_type('cell') %             -- Cell array
  spm_type('struct') %           -- Structure array
  spm_type('int8') %             -- 8-bit signed integer array
  spm_type('uint8') %            -- 8-bit unsigned integer array
  spm_type('int16') %            -- 16-bit signed integer array
  spm_type('uint16') %           -- 16-bit unsigned integer array
  spm_type('int32') %            -- 32-bit signed integer array
  spm_type('uint32') %           -- 32-bit unsigned integer array
  spm_type('int64') %            -- 64-bit signed integer array
  spm_type('uint64') %           -- 64-bit unsigned integer array
end

if isnan(spm_type(class(data)))
  % convert every unsupported data type into double
  data = double(data);
end

typ        = spm_type(class(data));
dim        = size(data);
Va         = [];
Va.mat     = transform;
Va.fname   = filename;
Va.dim     = [dim typ];
Va.n       = 1;
Va.pinfo   = [1 0 0]';
Va.private.hdr.dime.datatype = typ;
Va         = spm_create_vol(Va);
Va         = spm_write_vol(Va,data);

