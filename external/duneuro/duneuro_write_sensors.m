function duneuro_write_sensors(sens, filename)

if isfield(sens, 'elecpos')
  % treat as EEG for now, we should also consider sEEG and iEEG, which may need a different treatment
  write_duneuro_electrode_file(sens.elecpos, filename);
end

if isfield(sens, 'coilpos')
  % treat as meg
  assert(iscell(filename) && numel(filename)==2); % two filenames should be given, one for the coil positions, and one for the orientations
  write_duneuro_coil_file(sens.coilpos, filename{1});
  write_duneuro_projection_file(sens.coilori, filename{2});
end

function write_duneuro_coil_file(coil_loc, coil_filename)
% write_duneuro_coil_file(coil_loc, coil_filename)
% Write the electrode file for Duneuro application
% coil_loc : 3D Cartesien position of the coils (or integration point), Ncoil x 3

% Authors: Takfarinas MEDANI, December 2019;     

[filepath,name,ext] = fileparts(coil_filename);
if isempty(ext) || ~strcmp(ext,'.txt')
    ext = '.txt';
end
coil_filename = fullfile(filepath,[name,ext]);

fid = fopen(coil_filename, 'wt+');
fprintf(fid, '%d %d %d  \n', coil_loc');
fclose(fid); 

function write_duneuro_projection_file(coil_orientation,coilOrientation_filename)
%  write_duneuro_projection_file(coil_orientation,coilOrientation_filename)

% Authors: Takfarinas MEDANI, Decembre 2019;     

[filepath,name,ext] = fileparts(coilOrientation_filename);
if isempty(ext) || ~strcmp(ext,'.txt')
    ext = '.txt';
end
coilOrientation_filename = fullfile(filepath,[name,ext]);
% Nb_coils = size(coil_orientation, 1); 
% generate triedre orientation for each dipole
% coil_pos_orie = [kron(coil_pos,ones(3,1)), kron(ones(Nb_dipole,1), eye(3))];
% coil_orie = [ kron(ones(Nb_coils,1), eye(3))];
% coil_orie = repmat(([1 0 0 0 1 0 0 0 1]),Nb_coils,1);
fid = fopen(coilOrientation_filename, 'wt+');
% fprintf(fid, '%d %d %d %d %d %d \n', coil_pos_orie');
fprintf(fid, '%d %d %d  \n', coil_orientation');
fclose(fid); 

function write_duneuro_electrode_file(channel_loc, electrode_filename)
%write_duneuro_electrode_file(channel_loc, electrode_filename)
% Write the electrode file for Duneuro application
% channel_loc : 3D cartisien position of the electrodes, Nelec x 3

% Authors: Takfarinas MEDANI, August 2019;     

[filepath,name,ext] = fileparts(electrode_filename);
if isempty(ext) || ~strcmp(ext,'.txt')
    ext = '.txt';
end
electrode_filename = fullfile(filepath,[name,ext]);

fid = fopen(electrode_filename, 'wt+');
fprintf(fid, '%d %d %d  \n', channel_loc');
fclose(fid); 

