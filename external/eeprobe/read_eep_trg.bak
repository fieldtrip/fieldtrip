function [trg] = read_eep_trg(fn);

% READ_EEP_TRG reads triggers from an EEProbe *.trg file
%
% This function returns an Nx1 array with the N triggers
%
% [trg] = read_eep_trg(filename)
%
% trg(i).time   ... trigger latency in ms
% trg(i).offset ... byte offset
% trg(i).code   ... trigger code (string)
% trg(i).type   ... numeric value of trg.code
%
% where i as number between 1 and N (the number of triggers found)
%
% An EEProbe trigger file is formatted like
%   0.00195312 256
%   0.000     10350  __
%   17.033   2242926   1
%   20.535   2701934   5
%   21.096   2775406  13
%   21.098   2775662   8
%   ...
% where the first column is the trigger latency in seconds, the second
% column is the byte offset in the file and the third column is the triggercode. 
% The triggercode can be numeric or a string. The first line of the file contains the
% sample duration.
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_CNT, READ_EEP_REJ, READ_EEP_AVR
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

% $Log: not supported by cvs2svn $
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
% ANT Software BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

trg = [];

fid = fopen(fn, 'rb');
if fid<0
   return
end

header = fgetl(fid);
while ~feof(fid)
  tmp = fscanf(fid, '%f %d %s', 3);
  if ~isempty(tmp)
    new.time   = tmp(1)*1000;			% in ms
    new.offset = tmp(2)+1;			% offset 1
    new.code   = char(tmp(3:end));		% string
    new.type   = str2double(new.code);		% numeric
    trg = [trg; new];
  end
end

fclose(fid);  
