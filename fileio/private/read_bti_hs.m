function [output, firstIndexPoint] = read_hs_file( filename, outfile)

%read_hs_file Reads in BTI-Headshape files
%   filename: file with the headshape informations
%   outfile: if present, a ctf ".shape" file is written
%   output: if present, a 3xN matrix containing the headshape-points
%
%   (C) 2007 by Thomas Hartmann

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
