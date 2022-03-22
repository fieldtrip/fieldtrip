function json = ft_read_json(filename)

% FT_READ_JSON reads information from a JSON file and represents it as a MATLAB structure
%
% Use as
%   json = ft_read_json(filename)
%
% See also FT_WRITE_JSON, FT_READ_TSV, FT_WRITE_TSV, JSONDECODE, JSONENCODE

% Copyright (C) 2018-2022, Robert Oostenveld
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

ft_info('reading ''%s''\n', filename);

% Only check whether they are already present on the path
hasspm12up = ft_hastoolbox('SPM12UP');
hasjsonlab = ft_hastoolbox('jsonlab');

% The default is to use the native MATLAB implementation, which is fine for reading
% but not so optimal for writing.

if hasspm12up
  % use the SPM12 implementation
  json = spm_jsonread(filename);

elseif hasjsonlab
  % use the JSONLAB implementation
  json = loadjson(filename);

else
  % use the native MATLAB implementation
  fid = fopen_or_error(filename, 'rt');
  jsontxt = fread(fid, [1 inf], 'char=>char');
  fclose(fid);
  json = jsondecode(jsontxt);
end

% convert all strings into char-arrays
json = ft_struct2char(json);
