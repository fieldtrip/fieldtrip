function [isinside, pt, coord] = linextriangle(p0, p1, plane)
%  [isinside,pt,coord]=linextriangle(p0,p1,plane)
%
%  calculate the intersection of a 3d line (passing two points)
%  with a plane (determined by 3 points)
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      p0: a 3d point in form of (x,y,z)
%      p1: another 3d point in form of (x,y,z), p0 and p1 determins the line
%      plane: a 3x3 matrix, each row is a 3d point in form of (x,y,z)
%             this is used to define a plane
% outputs:
%      isinside: a boolean variable, 1 for the intersection is within the
%               3d triangle determined by the 3 points in plane; 0 is outside
%      pt: the coordinates of the intersection pint
%      coord: 1x3 vector, if isinside=1, coord will record the barycentric
%          coordinate for the intersection point within the triangle;
%          otherwise it will be all zeros.
%
% for degenerated lines or triangles, this will stop
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

[a, b, c, d] = getplanefrom3pt(plane);

if (a * a + b * b + c * c == 0.0)
    error('degenerated plane');
end

dl_n = sum([a b c] .* (p1 - p0));

if (dl_n == 0.0)
    error('degenerated line');
end

% solve for the intersection point
t = -(a * p0(1) + b * p0(2) + c * p0(3) + d) / dl_n;
pt = p0 + (p1 - p0) * t;

dist = sum(abs(diff(plane)));
[md, imax] = sort(dist);
if (md(2) == 0.0)
    error('degenerated triangle');
end
goodidx = imax(2:end);

ptproj = pt(goodidx);
mat0 = [plane(:, goodidx)', ptproj'; 1 1 1 1];

isinside = 0;
coord = [0 0 0];

det1 = det(mat0(:, [4 2 3], :));
det2 = det(mat0(:, [1 4 3], :));
if (det1 * det2 < 0)
    return
end
det3 = det(mat0(:, [1 2 4], :));
if (det2 * det3 < 0)
    return
end
if (det1 * det3 < 0)
    return
end
isinside = 1;
det0 = det(mat0(:, 1:3));

coord = [det1 det2 det3] / det0;
