function duneuro_write_sensors(sens, filename)

if isfield(sens, 'elecpos')
  % treat as EEG for now, we should also consider sEEG and iEEG, which may need a different treatment
  write_duneuro_nx3(sens.elecpos, filename);
end

if isfield(sens, 'coilpos')
  % treat as meg
  assert(iscell(filename) && numel(filename)==2); % two filenames should be given, one for the coil positions, and one for the orientations
  write_duneuro_nx3(sens.coilpos, filename{1});
  write_duneuro_nx3(sens.coilori, filename{2});
end

function write_duneuro_nx3(pos, filename)
% write_duneuro_nx3(coilpos, filename)
% Write a Nx3 text file for Duneuro application
% coilpos/coilori/elecpos: 3D Cartesian positions of the coils(or integration points)/coil orientations/electrodes, N x 3
% Authors: Takfarinas MEDANI, December 2019;

[filepath,name,ext] = fileparts(filename);
if isempty(ext) || ~strcmp(ext,'.txt')
  ext = '.txt';
end
filename = fullfile(filepath,[name,ext]);

fid = fopen(filename, 'wt+');
fprintf(fid, '%d %d %d  \n', pos');
fclose(fid);
