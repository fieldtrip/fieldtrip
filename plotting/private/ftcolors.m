function rgb = ftcolors(N)

% FTCOLORS returns an Nx3 rgb matrix with the
% colors of the fieldtrip logo at its extremes.
% Can be used as a colormap by FT_COLORMAP
%
% Use as:
%   rgb = ftcolors(N), or
%   rgb = ftcolors
%
% Without input arguments, N will be set to 64

if nargin==0
  N = 64;
end

top    = [172 27 42];
bottom = [18 140 120];
middle = [255 255 255];

if mod(N,2)==0
  N = N+1;
  dointerp = true;
end

n = (N+1)/2;
  
rgb = [linspace(bottom(1),middle(1),n) middle(1) linspace(middle(1),top(1),n); ...
       linspace(bottom(2),middle(2),n) middle(2) linspace(middle(2),top(2),n); ...
       linspace(bottom(3),middle(3),n) middle(3) linspace(middle(3),top(3),n)];
rgb(:,[n n+2]) = [];

if dointerp
  x = 1:N;
  xnew = 1.5:1:N-0.5;
  rgbnew(1,:) = interp1(x,rgb(1,:),xnew,'linear');
  rgbnew(2,:) = interp1(x,rgb(2,:),xnew,'linear');
  rgbnew(3,:) = interp1(x,rgb(3,:),xnew,'linear');
  rgb = rgbnew;
end

rgb = rgb'./255;

