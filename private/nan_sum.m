% nan_sum() - Take the sum, not considering NaN values
%
% Usage: same as sum()

% Author: Arnaud Delorme, CNL / Salk Institute, 16 Oct 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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


function out = nan_sum(in, dim)

    if nargin < 1
        help nan_sum;
        return;
    end;
    if nargin < 2
        if size(in,1) ~= 1
            dim = 1;
        elseif size(in,2) ~= 1
            dim = 2;
        else 
            dim = 3; 
        end;
    end;
    tmpin = in;
    tmpin(find(isnan(in(:)))) = 0;
    out = sum(tmpin, dim);
    
