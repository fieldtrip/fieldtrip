
function map = mapcolor2(A, breaks)
%MAPCOLOR2 Create a blue-gray-red colormap.
%   MAPCOLOR2(A, BREAKS), when A is a matrix and BREAKS a vector, returns a 
%   colormap that can be used as a parameter in COLORMAP function. BREAKS
%   has to contain six break values for a 7-class colormap. The break
%   values define the points at which the scheme changes color. The first 
%   three colors are blue, the middle one gray and the last three ones
%   are red.
%  
%   Example: A = ones(100, 100);
%            for i=1:10:100
%               A(i:end, i:end) = A(i,i) + 5;
%            end
%            A = A + rand(100, 100);
%            map = mapcolor2(A, [10 15 20 25 30 35]);  
%            pcolor(A), shading flat
%            colormap(map), colorbar
%

% Copyright (c) 2006 Markus Siivola

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% there must be 6 break values for seven colors
if length(breaks) ~= 6
   error('Break point vector must have 6 elements.');
end

% sort break values in ascending order if necessary
if ~issorted(breaks)
    breaks = sort(breaks);
end

% check the value range consumed by the data
rng = [min(A(:)) max(A(:))];

if any(breaks > rng(2)) | any(breaks < rng(1))
    error('Break point out of value range');
end

if ispc
    n = 256;
else
    n = 1000;
end

% a color scheme from blue through gray to red
colors = flipud([ 140 0   0;
                  194 80  68;
                  204 143 151;
                  191 191 191;
                  128 142 207;
                  70  84  158;
                  7   39  115 ] / 255);

% create a colormap with 256 colors
map = []; n = 256; ibeg = 1;
for i = 1:6
    iend = round(n * (breaks(i) - rng(1)) / (rng(2) - rng(1)));
    map  = [map; repmat(colors(i,:), [iend-ibeg+1 1])];
    ibeg = iend+1;
end
map = [map; repmat(colors(7,:), [n - size(map, 1) 1])];


