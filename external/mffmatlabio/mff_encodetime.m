% mff_encodetime - encode EGI time format from Matlab time format
%
% Usage:
%   timeout = mff_encodetime(timein);
%
% Input:
%  timein  - EGI time format (numerical)
%
% Output:
%  timeout - EGI string format time
%
% Example:
% timeValOut = mff_decodetime('2009-04-16T21:52:47.893250-08:00')
% mff_encodetime(timeValOut,'08:00')
%
% timeValOut = mff_decodetime('2009-04-16T21:52:47.893750-08:00')
% mff_encodetime(timeValOut,'08:00') % issue here

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function timeValOut = mff_encodetime(timeVal, timeZone)

    % remove GMT part
    microSec   = 0.000000000011525116860866546630859375000;
    
    % estimate error
    tmp = datestr(timeVal, 'yyyy-mm-ddTHH:MM:SS.FFF');
    remain = timeVal-mff_decodetime([tmp '000']);
    if remain < 0 
        tmp = datestr(timeVal-microSec*1000, 'yyyy-mm-ddTHH:MM:SS.FFF');
        remain = timeVal-mff_decodetime([tmp '000']);
        if remain < 0
            error('Negative microseconds');
        end
    end
    timeValOut = [ tmp sprintf('%.3d', round(remain/microSec/10)*10) timeZone ];

    
% old code
%     millisec   = 0.000000000011525116860866546630859375000;
%     
%     tmp = timeVal/millisec;
%     remain = (tmp-floor(tmp))*1000;
%     
%     timeValOut = [ datestr(timeVal, 'yyyy-mm-ddTHH:MM:SS.FFF') sprintf('%.3d', round(remain,-1)) '-' timeZone ];
