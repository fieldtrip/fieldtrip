% nan_mean() - Average, not considering NaN values
%
% Usage: same as mean()

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

% $Log: nan_mean.m,v $
% Revision 1.3  2006/03/14 15:01:58  roboos
% turn divideByZero warning temporary off
%
% Revision 1.2  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.1  2004/09/27 15:20:08  roboos
% subfunction required by multiplotTFR, singleplotTFR and topoplotTFR
%
% Revision 1.2  2002/10/17 18:43:16  arno
% debugging dim
%
% Revision 1.1  2002/10/17 02:34:52  arno
% Initial revision
%

function out = nan_mean(in, dim)

    if nargin < 1
        help nan_mean;
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
    warning('off', 'MATLAB:divideByZero');
    out = sum(tmpin, dim) ./ sum(~isnan(in),dim);
    warning('on', 'MATLAB:divideByZero');
    
