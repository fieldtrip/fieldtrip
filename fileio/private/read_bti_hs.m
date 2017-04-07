function [output, firstIndexPoint] = read_4d_hs( filename, outfile)

%read_hs_file Reads in BTI-Headshape files
%   filename: file with the headshape informations
%   outfile: if present, a ctf ".shape" file is written
%   output: if present, a 3xN matrix containing the headshape-points
%
%   (C) 2007 by Thomas Hartmann

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

if nargin == 1
    outfile = [];
end %if

fid       = fopen(filename, 'r', 'b');
version   = fread(fid, 1, '*uint32');
timestamp = fread(fid, 1, '*int32');
checksum  = fread(fid, 1, '*int32');
nPoints   = fread(fid, 1, '*int32');

firstIndexPoint = fread(fid, [3, 5], 'double')';

points = fread(fid, [3, double(nPoints)], 'double');

fclose(fid);

if(nargout > 0)
    output = points';
end %if

if(nargin == 2)
    fid = fopen(outfile, 'wt');
    fprintf(fid, '%d\n', nPoints);
    for i = 1:size(points, 2)
        fprintf(fid, '%.3f\t%.3f\t%.3f\n', points(1, i), points(2, i), points(3, i));
    end %for
    fclose(fid);

end %if
