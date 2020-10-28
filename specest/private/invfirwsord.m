% invfirwsord() - Estimate windowed sinc FIR filter transition band width
%                 depending on filter order and window type 
%
% Usage:
%   >> [df, dev] = invfirwsord(wtype, fs, m);
%   >> df = invfirwsord('kaiser', fs, m, dev);
%
% Inputs:
%   wtype - char array window type. 'rectangular', 'bartlett', 'hann',
%           'hamming', 'blackman', or 'kaiser'
%   fs    - scalar sampling frequency}
%   m     - scalar filter order
%   dev   - scalar maximum passband deviation/ripple (Kaiser window
%           only)
%
% Output:
%   df    - scalar estimated transition band width
%   dev   - scalar maximum passband deviation/ripple
%
% References:
%   [1] Smith, S. W. (1999). The scientist and engineer's guide to
%       digital signal processing (2nd ed.). San Diego, CA: California
%       Technical Publishing.
%   [2] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal
%       Processing: Principles, Algorithms, and Applications (3rd ed.).
%       Englewood Cliffs, NJ: Prentice-Hall
%   [3] Ifeachor E. C., & Jervis B. W. (1993). Digital Signal
%       Processing: A Practical Approach. Wokingham, UK: Addison-Wesley
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firws, firwsord

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005-2014 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Id$

function [ df, dev ] = invfirwsord(wintype, fs, m, dev)

winTypeArray = {'rectangular', 'bartlett', 'hann', 'hamming', 'blackman', 'kaiser'};
winDfArray = [0.9 2.9 3.1 3.3 5.5];
winDevArray = [0.089 0.056 0.0063 0.0022 0.0002];

% Check arguments
if nargin < 3 || isempty(fs) || isempty(m) || isempty(wintype)
    ft_error('Not enough input arguments.')
end

% Window type
wintype = find(strcmp(wintype, winTypeArray));
if isempty(wintype)
    ft_error('Unknown window type.')
end

if wintype == 6 % Kaiser window
    if nargin < 4 || isempty(dev)
        ft_error('Not enough input arguments.')
    end
    devdb = -20 * log10(dev);
    df = (devdb - 8) / (2.285 * 2 * pi * (m - 1));
else
    df = winDfArray(wintype) / m;
    dev = winDevArray(wintype);
end

% df is normalized
df = df * fs;

end
