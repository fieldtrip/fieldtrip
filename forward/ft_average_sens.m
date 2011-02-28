function [asens, afiducials] = ft_average_sens(sens, varargin)

% FT_AVERAGE_SENS computes average sensor array from a series of input
% arrays. Corresponding average fiducials can also be computed (optional)
%
% Use as
%   [asens, afid] = ft_average_sens(sens)
% where sens is a 1xN structure array containing N sensor arrays
%
% Additional options should be specified in key-value pairs and can be
%   'weights'    a vector of weights (will be normalized to sum==1)
%   'fiducials'  optional structure array of headshapes
%
% See also FT_READ_SENS, FT_PREPARE_VOL_SENS, FT_TRANSFORM_SENS

% Copyright (C) 2008-2011, Robert Oostenveld & Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% get the options
% fileformat = keyval('fileformat',  varargin);
weights       = keyval('weights',   varargin);
fiducials     = keyval('fiducials', varargin);
fb            = keyval('feedback',  varargin);

toplot     = ~isempty(fb);
nsens      = numel(sens);

% set the defaults
if isempty(weights) || ~any(weights)
  weights = ones(1, nsens);
end

% ensure that they are positive and sum to one
weights = abs(weights)./sum(abs(weights));

% determine whether the input contains EEG or MEG sensors
iseeg = ft_senstype(sens(1), 'eeg');
ismeg = ft_senstype(sens(1), 'meg');
istra = isfield(sens, 'tra');

if nsens == 1
  % no need to do any averaging and transformation
  asens      = sens;
  afiducials = fiducials;
  return;
end

pnt = detrend(sens(1).pnt, 'constant');
[u s v]= svd(pnt'*pnt);

% starting with the second PC works a little better (e.g. on BTi)
x1 = pnt*u(:, 2);
x2 = pnt*u(:, 1);

% detemine the indices of three reference sensors that are close to the three principal axes
[m ind1] = min(x1);
[m ind2] = max(x1);
[m ind3] = max(abs(x2 - mean(x2([ind1 ind2]))));

% compute the mean of the three reference sensors for the principal axes
mean1 = [0 0 0];
mean2 = [0 0 0];
mean3 = [0 0 0];

if istra
  tra = zeros(size(sens(1).tra));
end

if toplot
  axes(fb);
end

for i=1:nsens
  if ~isequal(sens(i).label, sens(1).label)
    error('all sensor arrays should have the same sensors for averaging');
  end
  
  % add them to the sum, to compute the mean location of each ref sensor
  mean1 = mean1 + weights(i).*sens(i).pnt(ind1, :);
  mean2 = mean2 + weights(i).*sens(i).pnt(ind2, :);
  mean3 = mean3 + weights(i).*sens(i).pnt(ind3, :);
  
  if istra
    tra = tra + weights(i).*sens(i).tra;
  end
  
  if toplot && ismeg
    plot3(sens(i).pnt(:, 1), sens(i).pnt(:, 2), sens(i).pnt(:, 3), '.', 'Color', [0.5 0.5 0.5]);
    hold on
  end
end

% compute the homogenous transformation matrix that aligns the three reference sensors to the x-, y-, and z-axis
tra = headcoordinates(mean1, mean2, mean3);

if ismeg
  % just realign the MEG coils
  tra1  = headcoordinates(sens(1).pnt(ind1, :), sens(1).pnt(ind2, :), sens(1).pnt(ind3, :));
  asens = ft_transform_sens(tra\tra1, sens(1));
  
elseif iseeg
  % also average sensor locations
  pnt = zeros(size(sens(1).pnt));
  for i=1:nsens
    tra1  = headcoordinates(sens(i).pnt(ind1, :), sens(i).pnt(ind2, :), sens(i).pnt(ind3, :));
    csens = ft_transform_sens(tra\tra1, sens(i));
    
    if toplot && iseeg
      plot3(csens.pnt(:, 1), csens.pnt(:, 2), csens.pnt(:, 3), '.', 'Color', [0.5 0.5 0.5]);
      hold on
    end
    
    pnt = pnt + weights(i).*csens.pnt;
  end % for nsens
  
  csens.pnt = pnt;
  asens     = csens;
  
else
  error('unsupported sensor type');
end

if istra
  % include the weighted averaged linear matrix that combines coils/electrodes into channels
  asens.tra = tra;
end

if toplot
  plot3(asens.pnt(:, 1), asens.pnt(:, 2), asens.pnt(:, 3), 'r.', 'MarkerSize', 20);
  hold on
  axis equal
  axis off
end

% average the fiducials
nfid = numel(fiducials);
switch nfid
  case 0
    % do nothing
    afiducials = [];
    
  case 1
    % FIXME this behaves inconsistent for eeg sensors, because tra1 happens to be the transformation for the last electrode set
    afiducials = ft_transform_headshape(tra\tra1, fiducials);
    
  case nsens
    fidpnt = zeros(size(fiducials(1).fid.pnt));
    
    for i=1:nsens
      if ~isequal(fiducials(i).fid.label, fiducials(1).fid.label)
        error('all fiducials should have the same labels for averaging');
      end
      
      tra1 = headcoordinates(sens(i).pnt(ind1, :), sens(i).pnt(ind2, :), sens(i).pnt(ind3, :));
      
      cfiducials = ft_transform_headshape(tra1, fiducials(i));
      
      fidpnt = fidpnt + weights(i).*cfiducials.fid.pnt;
    end
    
    cfiducials.pnt = [];
    cfiducials.fid.pnt = fidpnt;
    afiducials = ft_transform_headshape(inv(tra), cfiducials);
    
  otherwise
    error('there should be either one set of fiducials or a set per sensor array');
end % switch nfid
