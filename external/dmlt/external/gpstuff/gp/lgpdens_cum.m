function [cum] = lgpdens_cum(bb,x1,x2)
% [CUM] = LGPDENS_CUM(BB)
%
% Description
%   Given Bayesian Bootstrap estimated density, integrates it from point x1
%   to point x2

% Copyright (c) 2012 Ernesto Ulloa
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.addRequired('bb',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('x1',@(x) ~isempty(x) && isreal(x))
ip.addRequired('x2', @(x) ~isempty(x) && isreal(x))
ip.parse(bb,x1,x2)

 [p,pq,xt]=lgpdens(bb);
 I1=min(find(xt>x1));
 I2=max(find(xt<x2));
 sd=xt(2)-xt(1);
 cum=sd*trapz(p(I1:I2));
 
end

