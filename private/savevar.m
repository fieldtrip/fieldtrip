function savevar(filename, varname, value, hashfile)

% SAVEVAR is a helper function for cfg.outputfile
%
% See also LOADVAR

% Copyright (C) 2010, Robert Oostenveld
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

assert(ischar(filename), 'file name should be a string');
assert(ischar(varname), 'variable name should be a string');

ft_info('writing ''%s'' to file ''%s''\n', varname, filename);

eval(sprintf('%s = value;', varname));

% Create path if it doesn't exist
if ~isfolder(fileparts(filename))
  mkdir(fileparts(filename))
end

s = whos(varname);

% if variable < ~500 MB, store it in old (uncompressed) format, which is faster
if (s.bytes < 500000000) && ~ft_platform_supports('matlabversion', -inf, '2019a')
  save(filename, varname, '-v7', '-nocompression');
else
  save(filename, varname, '-v7.3');
end

% Also store the hash of the data we just stored into the specified
% hashfile, if requested. This is used by the reproducescript functionality
% in order to match up different input and output variables among different
% FieldTrip function calls.
if nargin > 3 && ~isempty(hashfile)
  % load the hashes already present in the file, if it exists
  if exist(hashfile, 'file')
    hashes = load(hashfile);
  else
    hashes = struct();
  end
  % key for the hash is the filename (no path) with 'f' prepended, since variable
  % names cannot begin with a number
  [p, f, x] = fileparts(filename);
  hashkey = ['f' f]; % only use the file name
  hashes.(hashkey) = ft_hash(value);
  % store it in the hashfile
  save(hashfile, '-struct', 'hashes');
  ft_info('writing data hash for ''%s'' to file ''%s''\n', varname, hashfile);
end
