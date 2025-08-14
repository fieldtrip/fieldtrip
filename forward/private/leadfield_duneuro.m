function [lf] = leadfield_duneuro(pos, headmodel, sens, method)

% LEADFIELD_DUNEURO computes EEG/MEG leadfields for a set of given dipoles
% using the finite element method (FEM)
%
% [lf] = leadfield_duneuro(pos, headmodel);
%
% with input arguments
%   pos     a matrix of dipole positions
%           (there can be 'deep electrodes', too)
%   vol     contains a FE volume conductor (output of ft_prepare_vol_sens)
%   method  string defining the modality ('eeg' or 'meg)
% The output lf is the leadfield matrix of dimensions m (rows) x n*3 (columns)

cfg = [];
cfg.post_process                   = bool2str(headmodel.duneuro.post_process);
cfg.subtract_mean                  = bool2str(headmodel.duneuro.subtract_mean);
cfg.source_model.type              = headmodel.duneuro.source_model.type;
cfg.source_model.initialization    = headmodel.duneuro.source_model.initialization;
cfg.source_model.intorderadd       = num2str(headmodel.duneuro.source_model.intorderadd);
cfg.source_model.intorderadd_lb    = num2str(headmodel.duneuro.source_model.intorderadd_lb);
cfg.source_model.numberOfMoments   = num2str(headmodel.duneuro.source_model.numberOfMoments);
cfg.source_model.referenceLength   = num2str(headmodel.duneuro.source_model.referenceLength);
cfg.source_model.relaxationFactor  = num2str(headmodel.duneuro.source_model.relaxationFactor);
cfg.source_model.weightingExponent = num2str(headmodel.duneuro.source_model.weightingExponent);
cfg.source_model.restrict          = bool2str(headmodel.duneuro.source_model.restrict);
cfg.source_model.mixedMoments      = bool2str(headmodel.duneuro.source_model.mixedMoments);

index = repmat(1:size(pos,1),3,1);
index = index(:);
dipoles = [pos(index,:)'; repmat(eye(3),1,size(pos,1))];


if isfield(headmodel, 'driver')
  try
    % compute the leadfield matrix
    lf = headmodel.driver.(sprintf('apply_%s_transfer', method))(headmodel.(sprintf('%s_transfer', method)), dipoles, cfg);
  catch
    ft_warning('An error occurred while computing the leadfield with duneuro.');
    rethrow(lasterror);
  end
else
  % save the required dipole positions in a file
  filename = fullfile(headmodel.duneuro.outputpath,'dipoles.txt');
  fid = fopen(filename, 'wt+');
  fprintf(fid, '%d %d %d %d %d %d \n', dipoles);
  fclose(fid);

  headmodel.duneuro.filename_dipoles = filename;
  headmodel.duneuro.modality         = method;

  % write the configuration file for the application
  duneuro_write_minifile(headmodel.duneuro, headmodel.duneuro.minifile_filename);

  % system call
  system([headmodel.duneuro.application ' ' headmodel.duneuro.minifile_filename]);

  % load the leadfield
  headmodel.duneuro = duneuro_read_leadfield(headmodel.duneuro);
end

switch method
  case 'eeg'
    if isfield(headmodel, 'driver')
      try
        %compute lead field matrix
        lf = headmodel.driver.apply_eeg_transfer(headmodel.eeg_transfer, dipoles, cfg);
      catch
        warning('An error occurred while computing the leadfield with duneuro.');
        rethrow(lasterror)
      end
    else
      % save the required dipole positions in a file
      filename = fullfile(headmodel.duneuro.outputpath,'dipoles.txt');
      fid = fopen(filename, 'wt+');
      fprintf(fid, '%d %d %d %d %d %d \n', dipoles);
      fclose(fid);

      headmodel.duneuro.filename_dipoles = filename;
      headmodel.duneuro.modality         = 'eeg';
      
      % write the configuration file for the application
      duneuro_write_minifile(headmodel.duneuro, headmodel.duneuro.minifile_filename);

      % system call
      system([headmodel.duneuro.application ' ' headmodel.duneuro.minifile_filename]);

      % load the leadfield
      headmodel.duneuro = duneuro_read_leadfield(headmodel.duneuro);
      
      % post processing is done outside
      lf = headmodel.duneuro.eeg.lf;
    end

  case 'meg'
    if isfield(headmodel,  'driver')
      try
        % compute lead field matrix
        Bs = headmodel.driver.apply_meg_transfer(headmodel.meg_transfer, dipoles, cfg);
      catch
        warning('An error occurred while computing the leadfield with duneuro.');
        rethrow(lasterror)
      end
    else
      % save the required dipole positions in a file
      filename = fullfile(headmodel.duneuro.outputpath,'dipoles.txt');
      fid = fopen(filename, 'wt+');
      fprintf(fid, '%d %d %d %d %d %d \n', dipoles);
      fclose(fid);

      headmodel.duneuro.filename_dipoles = filename;
      headmodel.duneuro.modality         = 'meg';
      
      % write the configuration file for the application
      duneuro_write_minifile(headmodel.duneuro, headmodel.duneuro.minifile_filename);

      % system call
      system([headmodel.duneuro.application ' ' headmodel.duneuro.minifile_filename]);

      % load the leadfield
      headmodel.duneuro = duneuro_read_leadfield(headmodel.duneuro);
      
      % post processing is done outside the if clause
      Bs = headmodel.duneuro.meg.Bs;
    end
    
    % compute primary B-field analytically
    Bp = compute_B_primary(sens.coilpos, dipoles', sens.coilori);
    
    % permeability constant mu in si units
    mu = 4*pi*1e-7; %unit: Tm/A
      
    % compute full B-field
    lf = mu/(4*pi) * (Bp - Bs);
end

function output = bool2str(input)

if input
  output = 'true';
else
  output = 'false';
end


function [Bp] = compute_B_primary(coils, dipoles, projections)

% compute primary magnetic B-field analytically
%
% input:
% coils (Nx3 matrix)
% dipoles (Mx6 matrix)
% projections (Nx3) matrix)

% check input
if size(coils,2)~=3
  error('Column size of coils must be 3.')
end

if size(dipoles,2)~=6
  error('Column size of dipoles must be 6.')
end

if size(projections,2)~=3
  error('Column size of projections must be 3.')
end

% apply formula of Sarvas

dip_pos = dipoles(:,1:3);
dip_mom = dipoles(:,4:6);
Bp = zeros(size(coils,1), size(dipoles,1));
for i = 1:size(coils,1)
  for j = 1 : size(dip_pos,1)
    R = coils(i,:);
    R_0 = dip_pos(j,:);
    A = R - R_0;
    a = norm(A);
    aa = A./(a^3);
    
    BpHelp = cross(dip_mom(j,:),aa);
    Bp(i,j) = BpHelp * projections(i, :)'; % projection of the primary B-field along the coil orientations
  end
end
