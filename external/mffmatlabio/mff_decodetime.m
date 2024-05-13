% mff_decodetime - decode EGI time format into Matlab time format
%
% Usage:
%   timeout = mff_decodetime(timein);
%   timeout = mff_decodetime(timein, fullResFlag);
%
% Input:
%  timein  - EGI time format (string)
%  fullResFlag - [0|1] flag of 1 do not subtract 730000 (looses some
%                resolution). Default is 0 and 730000 is subtracted.
%
% Output:
%  timeout - Matlab numerical time (see datenum). Note that 730000 is
%            subtracted to achive 

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

function timeValOut = mff_decodetime(timeVal, fullResFlag)

    % remove GMT part
    microSec     = 0.0000000000115251168608665466308593750000000;
    if nargin < 2 || isequal(fullResFlag, 0)
        timeVal(1)   = '0'; % fix the precision issue
    end
    timeValOut   = datenum(timeVal(1:23), 'yyyy-mm-ddTHH:MM:SS.FFF') + str2double(timeVal(24:26))*microSec;

%     indDash = find(timeVal == '-');
%     if indDash(end) > 10 else indDash = length(timeVal)+1; end;
%     
%     microSeconds = timeVal(indDash(end)-4:indDash(end)-1);
%     microSeconds = str2double(microSeconds)/24/60/60/1000;
%     fprintf('%1.20f\n', timeValOut);
%     timeValOut   = timeValOut+microSeconds;
%     fprintf('%1.20f\n', timeValOut);
    
