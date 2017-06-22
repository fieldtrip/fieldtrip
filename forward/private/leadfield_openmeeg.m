function [lp, voxels_in] = leadfield_openmeeg ( voxels, vol, sens, varargin )

% FT_OM_COMPUTE_LEAD uses OpenMEEG to compute the lead fields / potentials
% using the boundary element method (BEM).
% The inputs are as follows:
% voxels        = an [Nx3] array of voxel locations.
% vol           = the volume structure containing bnd and cond fields. In 
%                 order to save the matrices computed by OpenMEEG, the
%                 fields 'path' and 'basefile' should also be provided.
%                 Matrices will be stored under the directory specified by
%                 'path' and 'basefile' will be used to generate the
%                 filename. If FT_OM_COMPUTE_LEAD is run with the same
%                 'path' and 'basefile' parameters and detects the
%                 corresponding files from the OpenMEEG process, it will
%                 use those files for further processing rather than
%                 creating/calculating them again.
% sens          = the sens structure (elec for EEG or grad for MEG)
%
% Additional parameters can be specified by setting the following vol
% fields:
% vol.ecog      = "yes"/["no"] allows the computation of lead potentials
%                 for ECoG grids. (sens must be an EEG elec structure)
% vol.method    = ["hminv"]/"adjoint" to use the adjoint method instead of 
%                 the inverse head matrix.
%
% By default, the BEM will be computed using the inverse head matrix
% method. This is slower than the adjoint method, but more efficient if the
% BEM needs to be computed multiple times when sensor positions have
% moved relative to fixed head coordinates. To implement this type of
% computation:
% 1) vol.path and vol.basefile should be specified so that OpenMEEG
%    matrices will be saved.
% 2) Once the first BEM is computed, copy the following files to a second
%    working directory:
%       - *_hm.bin
%       - *_hminv.bin
%       - *_dsm#.bin (where # is a number)
% 3) vol.path should be changed to the second working directory and
%    vol.basefile should remain the same.
% 4) Recompute the BEM with the new set of sensor positions and/or voxels.
%
% Copyright (C) 2013, Daniel D.E. Wong, Sarang S. Dalal
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
%

vol = fixpos(vol); % renames old subfield 'pnt' to 'pos', if necessary

% Variable declarations
CPU_LIM = feature('numCores');
VOXCHUNKSIZE = 30000; % if OpenMEEG uses too much memory for a given computer, try reducing VOXCHUNKSIZE
om_format = 'binary'; % note that OpenMEEG's mat-file supported is limited in file size (2GB?)
switch(om_format)
    case 'matlab'
        om_ext = '.mat';
    case 'binary'
        om_ext = '.bin';
    otherwise
        error('invalid OpenMEEG output type requested');
end

OPENMEEG_PATH = []; % '/usr/local/bin/';    % In case OpenMEEG executables omitted from PATH variable

persistent ldLibraryPath0;

if ispc
    warning('Sorry, Windows is not yet tested');
elseif isunix
    setenv('OMP_NUM_THREADS',num2str(CPU_LIM));
    if(~ismac) % MacOS doesn't use LD_LIBRARY_PATH; in case of problems, look into "DYLD_LIBRARY_PATH"
        if isempty(ldLibraryPath0)
            ldLibraryPath0 = getenv('LD_LIBRARY_PATH');          % We'll restore this at the end
        end
        UNIX_LDLIBRARYPATH = '/usr/lib:/usr/local/lib';
        setenv('LD_LIBRARY_PATH',UNIX_LDLIBRARYPATH);   % MATLAB changes the default LD_LIBRARY_PATH variable
    end
end

[om_status,om_errmsg] = system(fullfile(OPENMEEG_PATH,'om_assemble')); % returns 0 if om_assemble is not happy
if(om_status ~= 0)
    error([om_errmsg 'Unable to properly execute OpenMEEG. Please configure variable declarations and paths in this file as needed.']);
else
    clear om_status
end

% Extra options
method = ft_getopt(vol,'method','hminv');
ecog = ft_getopt(vol,'ecog','no');

% Use basefile and basepath for saving files
if isfield(vol,'basefile')
    basefile = vol.basefile;
else
    basefile = tempname;
end
if isfield(vol,'path')
    path = vol.path;
    cleanup_flag = false;
else
    path = fullfile(tempdir,'ft-om');
    cleanup_flag = true;
end
mkdir(path);

sensorFile = fullfile(path, [basefile '_sensorcoords.txt']);
if exist(sensorFile,'file')
    disp('Sensor coordinate file already exists. Skipping...')
else
    disp('Writing sensor coordinates...')
    fid = fopen(sensorFile,'w');
    if ft_senstype(sens, 'eeg')
        for ii=1:size(sens.chanpos,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\n', sens.label{ii}, sens.chanpos(ii,:));
        end
    else    % MEG
        % Find channel labels for each coil -- non-trivial for MEG gradiometers!
        % Note that each coil in a gradiometer pair will receive the same label
        [chanlabel_idx,coilpos_idx]=find(abs(sens.tra)==1);

        newchanlabelmethod = true
        if(newchanlabelmethod)
            for ii=1:size(sens.coilpos,1)
                fprintf(fid,'%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',sens.coilpos(ii,:),sens.coilori(ii,:));
            end
        else
            if(size(sens.tra,1) < max(chanlabel_idx) || size(sens.tra,2) ~= length(coilpos_idx) || length(coilpos_idx) ~= size(sens.coilpos,1))
                % These dimensions should match; if not, some channels may have been
                % removed, or there's unexpected handling of MEG reference coils
                error('Mismatch between number of rows in sens.tra and number of channels... possibly some channels removed or unexpected MEG reference coil configuration');
            end
            
            for ii=1:length(coilpos_idx)
                coilpair_idx = find(chanlabel_idx(ii) == chanlabel_idx);
                
                if(length(coilpair_idx)==2)
                    whichcoil = find(ii == coilpair_idx);
                    
                    switch(whichcoil)
                        case 1
                            labelsuffix = 'A';
                        case 2
                            labelsuffix = 'B';
                    end
                else
                    labelsuffix = '';
                end
                label = [sens.label{chanlabel_idx(ii)} labelsuffix];
                fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',label,sens.coilpos(ii,:),sens.coilori(ii,:));
            end
        end
    end
    fclose(fid);
end

condFile = fullfile(path, [basefile '.cond']);
if exist(condFile,'file')
    disp('Conductivity file already exists. Skipping...')
else
    disp('Writing conductivity file...')
    write_cond(vol,condFile);
end

geomFile = fullfile(path, [basefile '.geom']);
if exist(geomFile,'file')
    disp('Geometry descriptor file already exists. Skipping...')
else
    disp('Writing geometry descriptor file...')
    write_geom(vol,geomFile,basefile);
end

disp('Writing OpenMEEG mesh files...')

write_mesh(vol,path,basefile);

disp('Validating mesh...')
[om_status om_msg] = system([fullfile(OPENMEEG_PATH, 'om_check_geom'), ' -g ', geomFile])
if(om_status ~= 0) % status = 0 if successful
    error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
end

disp('Writing dipole file...')
chunks = ceil(size(voxels,1)/VOXCHUNKSIZE);
dipFile = cell(chunks,1);
for ii = 1:chunks
    dipFile{ii} = fullfile(path, [basefile '_voxels' num2str(ii) om_ext]);
    if exist(dipFile{ii},'file')
        fprintf('\t%s already exists. Skipping...\n', dipFile{ii});
    else
        voxidx = ((ii-1)*VOXCHUNKSIZE + 1) : (min((ii)*VOXCHUNKSIZE,size(voxels,1)));
        writevoxels = [kron(voxels(voxidx,:),ones(3,1)) , kron(ones(length(voxidx),1),eye(3))];
        om_save_full(writevoxels,dipFile{ii},om_format);
    end
end


hmFile = fullfile(path, [basefile '_hm' om_ext]);
if exist(hmFile,'file')
    disp('Head matrix already exists. Skipping...')
else
    disp('Building head matrix')
    [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -hm ', geomFile, ' ', condFile, ' ', hmFile])
    if(om_status ~= 0) % status = 0 if successful
        error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
    end
end

if strcmp(method,'hminv')
    hminvFile = fullfile(path, [basefile '_hminv' om_ext]);
    if exist(hminvFile,'file')
        disp('Inverse head matrix already exists. Skipping...');
    else
        disp('Computing inverse head matrix');
        if(CPU_LIM >= 4) % Matlab's inverse function is multithreaded and performs faster with at least 4 cores
            om_save_sym(inv(om_load_sym(hmFile,om_format)),hminvFile,om_format);
        else
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_minverser'), ' ', hmFile, ' ', hminvFile])
            if(om_status ~= 0) % status = 0 if successful
                error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
            end
        end
    end
end
    
dsmFile = cell(chunks,1);
for ii = 1:chunks
    dsmFile{ii} = fullfile(path, [basefile '_dsm' num2str(ii) om_ext]);
    if exist(dsmFile{ii},'file')
        fprintf('\t%s already exists. Skipping...\n', dsmFile{ii});
    else
        disp('Assembling source matrix');
        [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -dsm ', geomFile, ' ', condFile, ' ', dipFile{ii}, ' ' dsmFile{ii}])
        if(om_status ~= 0) % status = 0 if successful
            error([om_msg, 'Aborting OpenMEEG pipeline due to above error. If 4-layer BEM attempted, try 3-layer BEM (scalp, skull, brain).']);
        end
    end
end


disp('--------------------------------------')


if ft_senstype(sens, 'eeg')
    if strcmp(ecog,'yes')
        ohmicFile = fullfile(path, [basefile '_h2ecogm']);
        cmd = '-h2ecogm';
    else
        ohmicFile = fullfile(path, [basefile '_h2em' om_ext]);
        cmd = '-h2em';
    end
else
    ohmicFile = fullfile(path, [basefile '_h2mm' om_ext]);
    cmd = '-h2mm';
end
if exist(ohmicFile,'file')
    disp('Ohmic current file already exists. Skipping...')
else
    disp('Calculating Contribution of Ohmic Currents')
    [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' ', cmd, ' ', geomFile, ' ', condFile, ' ' , sensorFile, ' ' , ohmicFile])
    if(om_status ~= 0) % status = 0 if successful
        error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
    end
end

if ft_senstype(sens, 'meg')
    disp('Contribution of all sources to the MEG sensors')
    scFile = cell(chunks,1);
    for ii = 1:chunks
        scFile{ii} = fullfile(path, [basefile '_ds2mm' num2str(ii) om_ext]);
        if exist(scFile{ii},'file')
            fprintf('\t%s already exists. Skipping...\n',scFile{ii})
        else
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -ds2mm ', dipFile{ii} ,' ', sensorFile, ' ' , scFile{ii}])
            if(om_status ~= 0) % status = 0 if successful
                error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
            end
        end
    end
end

disp('Putting it all together.')
bemFile = cell(chunks,1);
for ii = 1:chunks
    if ft_senstype(sens, 'eeg')
        bemFile{ii} = fullfile(path, [basefile '_eeggain' num2str(ii) om_ext]);
    else
        bemFile{ii} = fullfile(path, [basefile '_meggain' num2str(ii) om_ext]);
    end
    
    if exist(bemFile{ii},'file')
        fprintf('/t%s already exists. Skipping...\n', bemFile{ii});
        continue;
    end
    
    if strcmp(method,'hminv')
        if ft_senstype(sens, 'eeg')
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -EEG ', hminvFile, ' ', dsmFile{ii}, ' ', ohmicFile, ' ', bemFile{ii}]);
        else
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -MEG ', hminvFile, ' ', dsmFile{ii}, ' ', ohmicFile,' ', scFile{ii}, ' ',bemFile{ii}]);
        end
    else    % Adjoint method
        if ft_senstype(sens, 'eeg')
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -EEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', bemFile{ii}]);
        else
            [om_status, om_msg] = system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -MEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', scFile{ii}, ' ',bemFile{ii}]);
        end
    end
    if(om_status ~= 0) % status = 0 if successful
        error([om_msg, 'Aborting OpenMEEG pipeline due to above error.']);
    end
end

% Import lead field/potential
[g, voxels_in] = import_gain(path, basefile, ft_senstype(sens, 'eeg'));
if (voxels_in ~= voxels) && (nargout == 1); warning('Imported voxels from OpenMEEG process not the same as function input.'); end;

lp = sens.tra*g; % Mchannels x (3 orientations x Nvoxels)

% Cleanup
if cleanup_flag
    rmdir(basepath,'s')
end
if (isunix && ~ismac)
    setenv('LD_LIBRARY_PATH',ldLibraryPath0);
end



function write_cond(vol,filename)
fid=fopen(filename,'w');
fprintf(fid,'# Properties Description 1.0 (Conductivities)\n');
tissues = {'Scalp\t%f\n'; 'Skull\t%f\n'; 'CSF\t%f\n'; 'Brain\t%f\n'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
fprintf(fid,'Air\t0\n');
for ii=1:length(vol.cond)
    fprintf(fid,tissues{ii},vol.cond(ii));
end
fclose(fid);


function write_geom(vol,filename,basepathfile)
fid=fopen(filename,'w');
fprintf(fid,'# Domain Description 1.0\n');
fprintf(fid,'Interfaces %i Mesh\n',length(vol.cond));
tissues={'_scalp.tri\n'; '_skull.tri\n'; '_csf.tri\n'; '_brain.tri\n'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
for ii = 1:length(vol.cond)
    fprintf(fid,[basepathfile tissues{ii}]);
end
fprintf('\n');
fprintf(fid,'Domains %i\n',length(vol.cond)+1);
domains={'Scalp'; 'Skull'; 'CSF'; 'Brain'};
if length(vol.cond)==3; domains = domains([1 2 4]); end;
fprintf(fid,'Domain Air %i\n',1);
for ii = 1:length(vol.cond)
    if ii < length(vol.cond)
        fprintf(fid,['Domain ' domains{ii} ' %i -%i\n'],ii+1,ii);
    else
        fprintf(fid,['Domain ' domains{ii} ' -%i\n'],ii);
    end
end
fclose(fid);


function write_mesh(vol,path,basefile)
tissues={'_scalp'; '_skull'; '_csf'; '_brain'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
for ii = 1:length(vol.cond)
    meshFile = fullfile(path, [basefile tissues{ii} '.tri']);
    if exist(meshFile,'file')
        fprintf('\t%s already exists. Skipping...\n', meshFile);
        break;
    else
        om_save_tri(fullfile(path, [basefile tissues{ii} '.tri']),vol.bnd(ii).pos,vol.bnd(ii).tri);
        %savemesh([path basefile tissues{ii} '.mesh'],vol.bnd(ii).vertices/1000,vol.bnd(ii).faces-1,-vol.bnd(ii).normals);
    end
end


function [g, voxels] = import_gain(path, basefile, eegflag)
om_format = 'binary';
switch(om_format)
    case 'matlab'
        om_ext = '.mat';
    case 'binary'
        om_ext = '.bin';
    otherwise
        error('invalid OpenMEEG output type requested');
end

if eegflag
    omgainfiles = dir(fullfile(path, [basefile '_eeggain*' om_ext]));
else
    omgainfiles = dir(fullfile(path, [basefile '_meggain*' om_ext]));
end
omvoxfiles = dir(fullfile(path, [basefile '_voxels*' om_ext]));

g=[];
voxels=[];

% join gain/voxel files
% [openmeeg calculation may have been split for memory reasons]
for ii=1:length(omgainfiles)
    g = [g om_load_full(fullfile(path, omgainfiles(ii).name),om_format)];
    voxels = [voxels;om_load_full(fullfile(path, omvoxfiles(ii).name),om_format)];
end
voxels = voxels(1:3:end,1:3);
