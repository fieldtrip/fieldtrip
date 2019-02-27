function [map rng] = mapcolor_log(A, x, s)
% MAPCOLOR(A, X) returns a colormap ranging from blue through gray
% to red. The middle gray range starts from X(1) and ends at X(2).
%
% In    A   =   data matrix (vector with min and max values is enough)
%       x   =   gray range
%       s   =   saturation points, s(1) for blues and s(2) for reds
%
% Out   map =   colormap of size 256x3
%       rng =   min and max values of input A
% 
% Example usage:
%       datavector contains values around 1, values between 0.8 and 1.2 are
%       set to gray and values larger than 3 are set to darkest red
%
%       map = mapcolor_saturation(data, [.8 1.2], [0 3]);
%       colormap(map);
%       colorbar;
% 
% Example with pcolor: 
%       A = 10*rand(10,10);
%       map = mapcolor(A, [3, 7]);
%       pcolor(A), shading flat
%       colormap(map), colorbar

% Copyright (c) 2006 Markus Siivola

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

rng = [min(A(:)) max(A(:))];
if nargin < 2
    m = mean(A(find(A > 0)));
    x = [m m]; 
end
if nargin < 3
    s = rng;
end

red1  = log10([140 1 1]/255); % dark red
% red2  = log10([194 80 68]/255); 
% red3  = log10([204 143 151]/255);
gray  = log10([191 191 191]/255); % light gray
% blue3 = log10([128 142 207]/255);
% blue2 = log10([70 84 158]/255);
blue1 = log10([7 39 115]/255); % dark blue

% Find the first and the last row in the colormap matrix
% that should contain gray color
n = 256;  % colors

i1 = round(n * (x(1) - rng(1)) / (rng(2) - rng(1)));
i2 = round(n * (x(2) - rng(1)) / (rng(2) - rng(1)));
grays = repmat(10.^gray, i2-i1+1, 1);

nbs = round((s(1)-rng(1))/(x(1)-rng(1))*(i1-1)); %saturated blue
nblue = i1-1-max(nbs,1); 
blues = [logspace(blue1(1), blue1(1), nbs)' logspace(blue1(2), blue1(2), nbs)' logspace(blue1(3), blue1(3), nbs)';
         logspace(blue1(1), gray(1),  nblue)' logspace(blue1(2), gray(2),  nblue)' logspace(blue1(3), gray(3),  nblue)'];
     

nrs = round((rng(2)-s(2))/(rng(2)-x(2))*(n-i2)); %saturated red
nred = n-i2-max(nrs,1);
reds  = [logspace(red1(1), red1(1), nrs)' logspace(red1(2), red1(2), nrs)' logspace(red1(3), red1(3), nrs)';
         logspace(red1(1), gray(1), nred)' logspace(red1(2), gray(2), nred)' logspace(red1(3), gray(3), nred)'];

% Create the color map
map = [blues; grays; flipud(reds)];
