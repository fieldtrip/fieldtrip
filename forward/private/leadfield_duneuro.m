function [lf] = leadfield_duneuro(pos, vol)

% LEADFIELD_DUNEURO computes EEG leadfields for a set of given dipoles
% using the finite element method (FEM)
%
% [lf] = leadfield_duneuro(pos, vol);
%
% with input arguments
%   pos     a matrix of dipole positions
%           (there can be 'deep electrodes', too)
%   vol     contains a FE volume conductor (output of ft_prepare_vol_sens)
%
% The output lf is the leadfield matrix of dimensions m (rows) x n*3 (columns)


% Copyright (C) 2017, Sophie Schrader
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
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.


try
  %% compute lead field matrix
  lf = zeros(size(3*pos,2),size(vol.transfer,2));
  cfg = [];
  cfg.post_process = vol.post_process;
  cfg.subtract_mean = vol.subtract_mean;
  for i=1:size(pos, 1)
    dipoles =  [repmat(pos(i,:),3,1) diag([1.0,1.0,1.0])]';
    for j=1:size(dipoles, 2)
      t = vol.driver.apply_eeg_transfer(vol.transfer, dipoles(:,j), cfg);
      lf(j,:) = t;
    end
  end
  lf = lf';
catch
  warning('an error occurred while computing leadfield with duneuro');
  rethrow(lasterror)
end
