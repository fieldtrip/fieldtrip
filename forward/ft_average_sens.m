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

% get the optional input arguments
% fileformat = ft_getopt(varargin, 'fileformat');
weights       = ft_getopt(varargin, 'weights');
fiducials     = ft_getopt(varargin, 'fiducials');
fb            = ft_getopt(varargin, 'feedback');

toplot = ~isempty(fb);
nsens  = numel(sens);

% ensure the correct representation
for i=1:nsens
  newsens(i) = ft_datatype_sens(sens(i));
end
sens = newsens;
clear newsens;

fiducials = fixpos(fiducials);

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
else
  asens = [];
  afiducials = [];
end

pos = detrend(sens(1).chanpos, 'constant');
[u, s, v]= svd(pos'*pos);

% starting with the second PC works a little better (e.g. on BTi)
x1 = pos*u(:, 2);
x2 = pos*u(:, 1);

% detemine the indices of three reference sensors that are close to the three principal axes
[dum, ind1] = min(x1);
[dum, ind2] = max(x1);
[dum, ind3] = max(abs(x2 - mean(x2([ind1 ind2]))));

% compute the mean of the three reference sensors for the principal axes
mean1 = [0 0 0];
mean2 = [0 0 0];
mean3 = [0 0 0];

if istra
  asens.tra = zeros(size(sens(1).tra));
end

if toplot
  axes(fb);
end

for i=1:nsens
  if ~isequal(sens(i).label(:), sens(1).label(:))
    ft_error('all sensor arrays should have the same sensors for averaging');
  end
  
  % add them to the sum, to compute the mean location of each ref sensor
  mean1 = mean1 + weights(i).*sens(i).chanpos(ind1, :);
  mean2 = mean2 + weights(i).*sens(i).chanpos(ind2, :);
  mean3 = mean3 + weights(i).*sens(i).chanpos(ind3, :);
  
  if istra
    % include the weighted averaged linear matrix that combines coils/electrodes into channels
    asens.tra = asens.tra + weights(i).*sens(i).tra;
  end
  
  if toplot && ismeg
    plot3(sens(i).chanpos(:, 1), sens(i).chanpos(:, 2), sens(i).chanpos(:, 3), '.', 'Color', [0.5 0.5 0.5]);
    hold on
  end
end

% compute the homogeneous transformation matrix that aligns the three reference sensors to the x-, y-, and z-axis
tra = ft_headcoordinates(mean1, mean2, mean3);

if ismeg
  % just realign the MEG coils
  tra1  = ft_headcoordinates(sens(1).chanpos(ind1, :), sens(1).chanpos(ind2, :), sens(1).chanpos(ind3, :));
  asens = ft_transform_sens(tra\tra1, sens(1));
  
elseif iseeg
  % also average sensor locations
  pos = zeros(size(sens(1).chanpos));
  for i=1:nsens
    tra1  = ft_headcoordinates(sens(i).chanpos(ind1, :), sens(i).chanpos(ind2, :), sens(i).chanpos(ind3, :));
    csens = ft_transform_sens(tra\tra1, sens(i));
    
    if toplot && iseeg
      plot3(csens.chanpos(:, 1), csens.chanpos(:, 2), csens.chanpos(:, 3), '.', 'Color', [0.5 0.5 0.5]);
      hold on
    end
    
    pos = pos + weights(i).*csens.chanpos;
  end % for nsens
  
  csens.chanpos = pos;
  asens     = csens;
  
else
  ft_error('unsupported sensor type');
end

if toplot
  plot3(asens.chanpos(:, 1), asens.chanpos(:, 2), asens.chanpos(:, 3), 'r.', 'MarkerSize', 20);
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
    tra1  = ft_headcoordinates(sens(1).chanpos(ind1, :), sens(1).chanpos(ind2, :), sens(1).chanpos(ind3, :));
    afiducials = ft_transform_headshape(tra\tra1, fiducials);
    
  case nsens    
    hspos = [];
      
    for i=1:nsens
      hs = strmatch('headshape', fiducials(i).fid.label);
      fiducials(i).fid.label(hs)  = [];      
      
      fiducials(i).pos = [fiducials(i).pos; fiducials(i).fid.pos(hs, :)];
      
      fiducials(i).fid.pos(hs, :) = [];      
      
      if i == 1
          fidpos = zeros(size(fiducials(1).fid.pos));
      end
      
      if ~isequal(fiducials(i).fid.label, fiducials(1).fid.label)
        ft_error('all fiducials should have the same labels for averaging');
      end
      
      tra1 = ft_headcoordinates(sens(i).chanpos(ind1, :), sens(i).chanpos(ind2, :), sens(i).chanpos(ind3, :));
      
      cfiducials = ft_transform_headshape(tra1, fiducials(i));
      
      fidpos = fidpos + weights(i).*cfiducials.fid.pos;
      
      hspos = [hspos; cfiducials.pos];
    end
    
    cfiducials.pos = hspos;
    cfiducials.fid.pos = fidpos;
    afiducials = ft_transform_headshape(inv(tra), cfiducials);
    
    afiducials = ft_determine_units(afiducials);
    
    % remove redundant headshape points (3 cm precision)
    tolerance = 3;
    switch afiducials.unit
        case 'mm'
            c = 0.1;
        case 'cm'
            c = 1;
        case 'dm'
            c = 10;
        case 'm'
            c = 100;
    end
    
    [dum, ind] = unique(round(c*afiducials.pos/tolerance), 'rows');
    
    afiducials.pos = afiducials.pos(ind, :);
    
  otherwise
    ft_error('there should be either one set of fiducials or a set per sensor array');
end % switch nfid
