function [output, firstIndexPoint] = read_hs_file( filename, outfile)

%read_hs_file Reads in BTI-Headshape files
%   filename: file with the headshape informations
%   outfile: if present, a ctf ".shape" file is written
%   output: if present, a 3xN matrix containing the headshape-points
%
%   (C) 2007 by Thomas Hartmann

% $Log: read_bti_hs.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2008/05/11 16:29:57  vlalit
% Changed the function to also output fiducials
%
% Revision 1.2  2007/08/06 09:23:57  roboos
% removed debug output, transposed output to Nx3
%
% Revision 1.1  2007/08/06 09:03:42  roboos
% version as obtained from Nathan Weisz on 1-Aug-2007
%

if nargin == 1
    outfile = [];
end %if

fid = fopen(filename, 'r', 'b');
version = fread(fid, 1, '*uint32');
timestamp = fread(fid, 1, '*int32');
checksum = fread(fid, 1, '*int32');
nPoints = fread(fid, 1, '*int32');

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
