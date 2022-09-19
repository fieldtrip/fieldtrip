function level = artifact_level(dat, metric, mval, sd, connectivity)

% This function is shared between FT_REJECTVISUAL and FT_BADCHANNEL
%
% Use as
%   level = artifact_level(dat, metric, mval, sd, connectivity)
% where
%   dat           = nchan*ntime, data of a single trial
%   metric        = string
%   mval          = mean value over all trials
%   sd            = standard deviation over all trials
%   connectivity  = nchan*nchan connectivity matrix
% and
%   level         = nchan*1 vector with values

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
% Copyright (C) 2006-2021, Robert Oostenveld
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

if nargin<3
  mval = [];
end
if nargin<4
  sd = [];
end
if nargin<5
  connectivity = [];
end

if startsWith(metric, 'neighb')
  % cfg.neighbours is used in the calling function to compute the NxN Boolean matrix
  assert(~isempty(connectivity), 'this requires cfg.neighbours, see FT_PREPARE_NEIGHBOURS');
end

if contains(metric, 'zvalue')
  % The mval and sd should be computed over all trials, they are used here for computing
  % the z-value for each trial.
  assert(~isempty(mval));
  assert(~isempty(sd));
end

nchan = size(dat, 1);

switch metric
  case 'var'
    level = nanstd(dat, [], 2).^2;
  case 'db'
    level = 10*log10(nanstd(dat, [], 2).^2);
  case 'std'
    level = nanstd(dat, [], 2);
  case 'min'
    level = nanmin(dat, [], 2);
  case 'max'
    level = nanmax(dat, [], 2);
  case 'maxabs'
    level = nanmax(abs(dat), [], 2);
  case 'range'
    level = nanmax(dat, [], 2) - nanmin(dat, [], 2);
  case 'kurtosis'
    level = kurtosis(dat, [], 2);
  case '1/var'
    level = 1./(nanstd(dat, [], 2).^2);
  case 'zvalue'
    level = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
  case 'maxzvalue'
    level = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
    
  case 'neighbexpvar'
    % this results in a Nx1 vector
    level = nan(nchan, 1);
    
    % compute for each channel the amount of variance that can be explained by its neighbours
    % values close to 1 are good
    for i=1:nchan
      nb = find(connectivity(i,:));
      if isempty(nb)
        level(i) = 0;
      else
        % Compute the residual using the GLM
        %   y = beta * x + residual
        % where y is a row vector, not a column vector as common with fMRI
        x = dat(nb, :);
        y = dat(i,:);
        % remove the mean and divide by the standard deviation
        x = ft_preproc_standardize(x);
        y = ft_preproc_standardize(y);
        if all(x(:) == 0) || all(y(:) == 0)
          level(i) = 0;
        else
          residual = y - y / x * x;
          level(i) = 1 - sum(residual.^2) / sum(y.^2);
        end
      end
    end
    
  case 'neighbcorr'
    % this results in a NxN matrix
    level = nan(nchan, nchan);
    
    % compute the correlation between each channel and each of its neighbours
    % values close to 1 are good
    for i=1:nchan
      nb = find(connectivity(i,:));
      for j=nb(:)'
        level(i,j) = corr(dat(i,:)', dat(j,:)');
      end
    end
    
  case 'neighbstdratio'
    % this results in a NxN matrix
    level = nan(nchan, nchan);
    
    % compute the standard deviation for each channel
    chanstd = nan(nchan,1);
    for i=1:nchan
      chanstd(i) = std(dat(i,:));
    end
    
    % compute the relative difference in standard deviation of each channel with that of each of its neighbours
    % values close to 0 are good
    for i=1:nchan
      nb = find(connectivity(i,:));
      for j=nb(:)'
        level(i,j) = abs(chanstd(i)-chanstd(j)) / chanstd(j);
      end
    end
    
  otherwise
    ft_error('unsupported method');
end