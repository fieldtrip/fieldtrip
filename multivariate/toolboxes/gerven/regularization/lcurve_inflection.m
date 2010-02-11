function [idx,fit,err] = lcurve_inflection(lcurve)
% LCURVE_INFLECTION finds the inflection point of an lcurve by using a
% piecewise linear fit using two segements. It returns the index (idx) of the
% lcurve which minimizes the squared error (err) between the lcurve and the
% piecewise linear fit (fit).
%
% SEE ALSO:
% spfit.m
%
% Copyright (c) 2009, Marcel van Gerven
%
% $Log: lcurve_inflection.m,v $
%

idx = 0;
fit = [];
p = [];
err = inf;

ll = length(lcurve);
for i=2:(ll-2)
  
%   xb = [1 i ll];
%   sp = spfit(1:ll,lcurve,xb,2);  
%   y = ppval(sp,1:length(lcurve));
%   
%   ss = norm(lcurve-y,2).^2;

  % split the data into two segments
  d1 = lcurve(1:i);
  d2 = lcurve((i+1):ll);

  b1 = polyfit(1:length(d1),d1,1);
  v1 = polyval(b1,1:length(d1));
 
  b2 = polyfit(1:length(d2),d2,1);
  v2 = polyval(b2,1:length(d2));
 
  ss = norm(d1 - v1,2).^2 + norm(d2 - v2,2).^2;

  if ss < err
    err = ss;
    idx = i;
    fit = [v1 v2];
  end
  
end