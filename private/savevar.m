function savevar(filename, varname, value, hashfile)

% SAVEVAR is a helper function for cfg.outputfile
%
% See also LOADVAR

% Copyright (C) 2010, Robert Oostenveld
%
% $Id$

assert(ischar(filename), 'file name should be a string');
assert(ischar(varname), 'variable name should be a string');

ft_info('writing ''%s'' to file ''%s''\n', varname, filename);

eval(sprintf('%s = value;', varname));

s = whos(varname);

% if variable < ~500 MB, store it in old (uncompressed) format, which is faster
if (s.bytes < 500000000)
  save(filename, varname, '-v6');
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