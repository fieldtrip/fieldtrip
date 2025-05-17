function [leftpt, leftcurve, rightpt, rightcurve] = slicesurf3(node, elem, p1, p2, p3, step, minangle)
%
% [leftpt,leftcurve,rightpt,rightcurve]=slicesurf3(node,elem,p1,p2,p3,step,minangle)
%
% Slice a closed surface by a plane and extract the landmark nodes along
% the intersection between p1 and p3, then output into 2 segments: between
% p2 to p1 (left half), and p2 to p3 (right half)
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: an N x 3 array defining the 3-D positions of the mesh
%    elem: an N x 3 interger array specifying the surface triangle indices;
%    p1: 3D position of the start on the curve-of-interest
%    p2: 3D position of the middle on the curve-of-interest
%    p3: 3D position of the end on the curve-of-interest
%    step: (optional) a percentage (0-100) specifying the spacing of the
%        output landmark nodes; step=20 means the landmarks on the left
%        curve are spaced as 20% of the total lengths of the left-half, and
%        those on the right-curve are spaced at 20% of the right-half,
%        starting from p2.
%    minangle: (optional) a positive minangle will ask this function to
%        call polylinesimplify to remove sharp turns on the curve.
%
% output:
%    leftpt: the equal-spaced landmark nodes on the left-half (p2-p1)
%            intersection curve; spacing between these nodes are
%            (step% * length of the curve between p2-p1)
%    leftcurve: all nodes on the left-half (p2-p1) intersection curve
%    rightpt: the equal-spaced landmark nodes on the right-half (p2-p3)
%            intersection curve; spacing between these nodes are
%            (step% * length of the curve between p2-p3)
%    rightcurve: all nodes on the left-half (p2-p1) intersection curve
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

fullcurve = slicesurf(node, elem, [p1; p2; p3]);
if (nargin >= 7 && minangle > 0)
    fullcurve = polylinesimplify(fullcurve, minangle);
end

[fulllen, fullcurve] = polylinelen(fullcurve, p1, p3, p2);

[leftlen,  leftcurve] = polylinelen(fullcurve, p2, p1);
if (nargin >= 6)
    [idx, weight, leftpt] = polylineinterp(leftlen, sum(leftlen) * (step:step:(100 - step * 0.5)) * 0.01, leftcurve);
else
    leftpt = leftcurve;
end

if (nargout > 2)
    [rightlen, rightcurve] = polylinelen(fullcurve, p2, p3);
    if (nargin >= 6)
        [idx, weight, rightpt] = polylineinterp(rightlen, sum(rightlen) * (step:step:(100 - step * 0.5)) * 0.01, rightcurve);
    else
        rightpt = rightcurve;
    end
end
