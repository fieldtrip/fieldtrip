function h = ft_hash(s)

% FT_HASH computes a MD5 hash from a MATLAB variable or structure
%
% It will first try a hashing algorithm implemented as a mex file.
% If that fails, it falls back to a slower one that is based on Java.

% both CalcMD5 and DataHash are from Mathworks file exchange
ft_hastoolbox('fileexchange', 1);

try
  % this one uses a mex file
  if isstruct(s) || iscell(s) || istable(s)
    h = CalcMD5(mxSerialize(s));
  else
    h = CalcMD5(s);
  end
catch
  % this one uses Java
  h = DataHash(s);
end
