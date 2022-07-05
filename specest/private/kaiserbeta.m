% kaiserbeta() - Estimate Kaiser window beta
%
% Usage:
%   >> beta = pop_kaiserbeta(dev);
%
% Inputs:
%   dev       - scalar maximum passband deviation/ripple
%
% Output:
%   beta      - scalar Kaiser window beta
%
% References:
%   [1] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal
%       Processing: Principles, Algorithms, and Applications (3rd ed.).
%       Englewood Cliffs, NJ: Prentice-Hall
%
% Author: Andreas Widmann, University of Leipzig, 2005

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

function [ beta ] = kaiserbeta(dev)

devdb = -20 * log10(dev);

if devdb > 50
    beta = 0.1102 * (devdb - 8.7);
elseif devdb >= 21
    beta = 0.5842 * (devdb - 21)^0.4 + 0.07886 * (devdb - 21);
else
    beta = 0;
end

end
