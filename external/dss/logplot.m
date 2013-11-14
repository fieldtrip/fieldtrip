function logplot(x_orig, logscale)
% Normal plot with layered logscale plot
%   logplot(x_orig, logscale)
%     x_orig    Plotted data
%     logscale  Optional 2-element vector representing plot range
%               for cropping logarithmic plot area

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if size(x_orig,1)==1
  x=x_orig';
else
  x=x_orig;
end

if size(x,1)==0; return; end
i=[1:size(x,1)];

x(x==0)=NaN;
y=log(x);
y(isinf(y))=NaN;

if nargin<2
  offset = min(x)-min(y);
  scale = (max(x)-min(x))/(max(y)-min(y));
else
  offset = logscale(1)-min(y);
  scale = (logscale(2)-logscale(1))/(max(y)-min(y));
end
plot(i, x, '-', i, (y-min(y))*scale+min(x), ':');
