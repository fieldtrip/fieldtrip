function [lf] = leadfield_duneuro(pos, vol, method)

% LEADFIELD_DUNEURO computes EEG/MEG leadfields for a set of given dipoles
% using the finite element method (FEM)
%
% [lf] = leadfield_duneuro(pos, vol);
%
% with input arguments
%   pos     a matrix of dipole positions
%           (there can be 'deep electrodes', too)
%   vol     contains a FE volume conductor (output of ft_prepare_vol_sens)
%   method  string defining the modality ('eeg' or 'meg)
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

cfg = [];
cfg.post_process                   = vol.post_process;
cfg.subtract_mean                  = vol.subtract_mean;
cfg.source_model.type              = vol.forward;
cfg.source_model.initialization    = vol.initialization;
cfg.source_model.intorderadd       = vol.intorderadd;
cfg.source_model.intorderadd_lb    = vol.intorderadd_lb;
cfg.source_model.numberOfMoments   = vol.numberOfMoments;
cfg.source_model.referenceLength   = vol.referenceLength;
cfg.source_model.relaxationFactor  = vol.relaxationFactor;
cfg.source_model.restrict          = vol.restrict;
cfg.source_model.weightingExponent = vol.weightingExponent;
cfg.source_model.mixedMoments      = vol.mixedMoments;

index = repmat(1:size(pos,1),3,1);
index = index(:);
dipoles = [pos(index,:)'; repmat(eye(3),1,size(pos,1))];

switch method
  case 'eeg'
    try
      %compute lead field matrix
      lf = vol.driver.apply_eeg_transfer(vol.eeg_transfer, dipoles, cfg);
    catch
      warning('An error occurred while computing the leadfield with duneuro.');
      rethrow(lasterror)
    end
    
  case 'meg'
    try
      % compute lead field matrix
      lf = vol.driver.apply_meg_transfer(vol.meg_transfer, dipoles, cfg);
    catch
      warning('An error occurred while computing the leadfield with duneuro.');
      rethrow(lasterror)
    end
end
