% MEX is a Matlab compatibility shim wrapper for the MEX function. It exists
% because Octave's mex() function returns a status code instead of raising
% an error when the compilation fails, unlike Matlab's mex(), which
% errors on failed compilation. This wrapper smooths over the difference,
% and makes mex() raise errors on failure in Octave, too.
%
%
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

function mex(varargin)
  real_mex_dir = [matlabroot, '/share/octave/' version '/m/miscellaneous'];
  %status = builtin('mex', varargin{:});
  unwind_protect
    addpath(real_mex_dir, '-begin');
    status = mex(varargin{:});
    if status != 0
      error('mex() invocation failed');
    endif
  unwind_protect_cleanup
    % Shuffle it back to the end of the path
    addpath(real_mex_dir, '-end');
  end_unwind_protect
endfunction
