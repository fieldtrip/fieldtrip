function ft_mex (varargin)

% FT_MEX is a compatibility shim wrapper for the MEX function. It exists
% because Octave's mex() function returns a status code instead of raising
% an error when the compilation fails, unlike Matlab's mex(), which
% errors on failed compilation. This wrapper smooths over the difference,
% and always raises an error on failure, regardless of whether you're 
% running in Matlab or Octave.
%
% The signature for FT_MEX is exactly the same as MEX for the system you
% are running it on. (All its arguments are passed directly on to MEX.)
% FT_MEX is a drop-in replacement for MEX, and should be used in all places
% in FieldTrip where MEX would be called.
%
% See also: MEX

% Copyright (C) 2019, Andrew Janke
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

persistent is_octave
if isempty(is_octave)
  v = ver;
  is_octave = ismember('Octave', {v.Name});
end

if is_octave
  status = mex(varargin{:});
  if status != 0
    ft_error('mex() invocation failed.');
  end
else
  mex(varargin{:});
end
