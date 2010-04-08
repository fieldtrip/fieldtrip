function [rej] = read_eep_rej(fn);

% READ_EEP_REJ reads rejection marks from an EEProbe *.rej file
%
% This function returns a Nx2 matrix with the begin and end latency
% of N rejection marks. The latency is in miliseconds.
%
% rej = read_eep_rej(filename)
%
% An EEProbe rejection file is formatted like
%   0.0000-0.3640
%   2.4373-3.5471
%   ... 
% where rejection begin and end are given in seconds. This function 
% converts the latency in miliseconds.
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_CNT, READ_EEP_TRG, READ_EEP_AVR
%

% Copyright (C) 2002, Robert Oostenveld
%                     Aalborg University, Denmark
%                     http://www.smi.auc.dk/~roberto/
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

% $Log: read_eep_rej.m,v $
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:17:02  jwiskerke
% Added m-files without binary code in maple distribution.
%
% Revision 1.2  2003/10/24 13:34:41  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% Revision 1.1.1.1  2003/03/11 15:24:51  roberto
% updated help and copyrights
% ANT Software BV, The Netherlands, www.ant-software.nl / info@ant-software.nl
%

rej = [];

fid = fopen(fn, 'rb');
if fid<0
   return 
end
while ~feof(fid)
  tmp = fscanf(fid, '%f-%f', 2);
  if ~isempty(tmp)
    rej = [rej; tmp'];
  end
end

% convert to ms
rej = 1000*rej;

fclose(fid);  