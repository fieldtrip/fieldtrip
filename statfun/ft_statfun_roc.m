function [s, cfg] = ft_statfun_roc(cfg, dat, design)

% FT_STATFUN_ROC computes the area under the curve (AUC) of the Receiver Operator
% Characteristic (ROC). This is a measure of the separability of the data observed in
% two conditions. The AUC can be used for statistical testing whether the two
% conditions can be distinguished on the basis of the data.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_roc'
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%
% The labels for the independent variable should be specified as the number 1 and 2.
%
% Note that this statfun performs a one sided test in which condition "1" is assumed
% to be larger than condition "2". This function does not compute an analytic
% probability of condition "1" being larger than condition "2", but can be used in a
% randomization test, including clustering.
%
% A low-level example with 10 channel-time-frequency points and 1000 observations per
% condition goes like this:
%   dat1 = randn(10,1000) + 1;
%   dat2 = randn(10,1000);
%   design = [1*ones(1,1000) 2*ones(1,1000)];
%   stat = ft_statfun_roc([], [dat1 dat2], design);
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the defaults
cfg.ivar = ft_getopt(cfg, 'ivar', 1);

% on-the-fly log transformation is not supported any more, please use FT_MATH instead
cfg = ft_checkconfig(cfg, 'forbidden', 'logtransform');
cfg = ft_checkconfig(cfg, 'forbidden', 'numbins');

% start with a quick test to see whether there appear to be NaNs
if any(isnan(dat(1,:)))
  % exclude trials that contain NaNs for all observed data points
  sel    = all(isnan(dat),1);
  dat    = dat(:,~sel);
  design = design(:,~sel);
end

% logical indexing is faster than using find(...)
selA = (design(cfg.ivar,:)==1);
selB = (design(cfg.ivar,:)==2);
% select the data in the two classes
datA = dat(:, selA);
datB = dat(:, selB);

nobs = size(dat,1);
na   = size(datA,2);
nb   = size(datB,2);
auc  = zeros(nobs, 1);

for k = 1:nobs
  % compute the area under the curve for each channel/time/frequency
  a = datA(k,:);
  b = datB(k,:);

  % to speed up the AUC, the critical value is determined by the actual
  % values in class B, which also ensures a regular sampling of the False Alarms
  b = sort(b);

  ca = zeros(nb+1,1);
  ib = zeros(nb+1,1);
  % cb = zeros(nb+1,1);
  % ia = zeros(nb+1,1);

  for i=1:nb
    % for the first approach below, the critval could also be choosen based on e.g. linspace(min,max,n)
    critval = b(i);

    % for each of the two distributions, determine the number of correct and incorrect assignments given the critical value
    % ca(i) = sum(a>=critval);
    % ib(i) = sum(b>=critval);
    % cb(i) = sum(b<critval);
    % ia(i) = sum(a<critval);

    % this is a much faster approach, which works due to using the sorted values in b as the critical values
    ca(i) = sum(a>=critval);   % correct assignments to class A
    ib(i) = nb-i+1;            % incorrect assignments to class B
  end

  % add the end point
  ca(end) = 0;
  ib(end) = 0;
  % cb(end) = nb;
  % ia(end) = na;

  hits = ca/na;
  fa   = ib/nb;

  % the numerical integration is faster if the points are sorted
  hits = fliplr(hits);
  fa   = fliplr(fa);

  if false
    % this part is optional and should only be used when exploring the data
    figure
    plot(fa, hits, '.-')
    xlabel('false positive');
    ylabel('true positive');
    title('ROC-curve');
  end

  % compute the area under the curve using numerical integration
  auc(k) = numint(fa, hits);

end

% return the area under the curve as the statistic of interest
s.stat = auc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMINT computes a numerical integral of a set of sampled points using
% linear interpolation. Alugh the algorithm works for irregularly sampled points
% along the x-axis, it will perform best for regularly sampled points
%
% Use as
%   z = numint(x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = numint(x, y)
if ~all(diff(x)>=0)
  % ensure that the points are sorted along the x-axis
  [x, i] = sort(x);
  y = y(i);
end
n = length(x);
z = 0;
for i=1:(n-1)
  x0 = x(i);
  y0 = y(i);
  dx = x(i+1)-x(i);
  dy = y(i+1)-y(i);
  z = z + (y0 * dx) + (dy*dx/2);
end
