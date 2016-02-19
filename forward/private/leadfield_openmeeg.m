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

% Variable declarations
CPU_LIM = feature('numCores');
VOXCHUNKSIZE = 25000;

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

if system(fullfile(OPENMEEG_PATH,'om_assemble'))
    error('Unable to properly execute OpenMEEG. Please configure variable declarations and paths in this file as needed.');
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
        for ii=1:size(sens.coilpos,1)
            if ii <= length(sens.label)
                label = sens.label{ii};
            else
                label = ['XLABEL' num2str(ii-length(sens.label))]; % Make up a label name
            end
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',label,sens.coilpos(ii,:),sens.coilori(ii,:));
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
[stat res] = system([fullfile(OPENMEEG_PATH, 'om_check_geom'), ' -g ', geomFile]);
if stat > 0
    error([res, '...quitting']);
end

disp('Writing dipole file...')
chunks = ceil(size(voxels,1)/VOXCHUNKSIZE);
dipFile = cell(chunks,1);
for ii = 1:chunks
    dipFile{ii} = fullfile(path, [basefile '_voxels' num2str(ii) '.bin']);
    if exist(dipFile{ii},'file')
        fprintf('\t%s already exists. Skipping...\n', dipFile{ii});
    else
        voxidx = ((ii-1)*VOXCHUNKSIZE + 1) : (min((ii)*VOXCHUNKSIZE,size(voxels,1)));
        writevoxels = [kron(voxels(voxidx,:),ones(3,1)) , kron(ones(length(voxidx),1),eye(3))];
        om_save_full(writevoxels,dipFile{ii},'binary');
    end
end


disp('--------------------------------------')

hmFile = fullfile(path, [basefile '_hm.bin']);
if exist(hmFile,'file')
    disp('Head matrix already exists. Skipping...')
else
    disp('Building head matrix')
    system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -hm ', geomFile, ' ', condFile, ' ', hmFile]);
end

if strcmp(method,'hminv')
    hminvFile = fullfile(path, [basefile '_hminv.bin']);
    if exist(hminvFile,'file')
        disp('Inverse head matrix already exists. Skipping...');
    else
        disp('Computing inverse head matrix');
        system([fullfile(OPENMEEG_PATH, 'om_minverser'), ' ', hmFile, ' ', hminvFile]);
    end
    
    dsmFile = cell(chunks,1);
    for ii = 1:chunks
        dsmFile{ii} = fullfile(path, [basefile '_dsm' num2str(ii) '.bin']);
        if exist(dsmFile{ii},'file')
            fprintf('\t%s already exists. Skipping...\n', dsmFile{ii});
        else
            disp('Assembling source matrix');
            system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -dsm ', geomFile, ' ', condFile, ' ', dipFile{ii}, ' ' dsmFile{ii}]);
        end
    end
end

if ft_senstype(sens, 'eeg')
    if strcmp(ecog,'yes')
        ohmicFile = fullfile(path, [basefile '_h2ecogm']);
        cmd = '-h2ecogm';
    else
        ohmicFile = fullfile(path, [basefile '_h2em.bin']);
        cmd = '-h2em';
    end
else
    ohmicFile = fullfile(path, [basefile '_h2mm.bin']);
    cmd = '-h2mm';
end
if exist(ohmicFile,'file')
    disp('Ohmic current file already exists. Skipping...')
else
    disp('Calculating Contribution of Ohmic Currents')
    system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' ', cmd, ' ', geomFile, ' ', condFile, ' ' , sensorFile, ' ' , ohmicFile]);
end

if ft_senstype(sens, 'meg')
    disp('Contribution of all sources to the MEG sensors')
    scFile = cell(chunks,1);
    for ii = 1:chunks
        scFile{ii} = fullfile(path, [basefile '_ds2mm' num2str(ii) '.bin']);
        if exist(scFile{ii},'file')
            fprintf('\t%s already exists. Skipping...\n',scFile{ii})
        else
            system([fullfile(OPENMEEG_PATH, 'om_assemble'), ' -ds2mm ', dipFile{ii} ,' ', sensorFile, ' ' , scFile{ii}]);
        end
    end
end

disp('Putting it all together. Creating the BEM might take a while')
bemFile = cell(chunks,1);
for ii = 1:chunks
    if ft_senstype(sens, 'eeg')
        bemFile{ii} = fullfile(path, [basefile '_eeggain' num2str(ii) '.bin']);
    else
        bemFile{ii} = fullfile(path, [basefile '_meggain' num2str(ii) '.bin']);
    end
    
    if exist(bemFile{ii},'file')
        fprintf('/t%s already exists. Skipping...\n', bemFile{ii});
        continue;
    end
    
    if strcmp(method,'hminv')
        if ft_senstype(sens, 'eeg')
            system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -EEG ', hminvFile, ' ', dsmFile{ii}, ' ', ohmicFile, ' ', bemFile{ii}]);
        else
            system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -MEG ', hminvFile, ' ', dsmFile{ii}, ' ', ohmicFile,' ', scFile{ii}, ' ',bemFile{ii}]);
        end
    else    % Adjoint method
        if ft_senstype(sens, 'eeg')
            system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -EEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', bemFile{ii}]);
        else
            system([fullfile(OPENMEEG_PATH, 'om_gain'), ' -MEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', scFile{ii}, ' ',bemFile{ii}]);
        end
    end
end

% Import lead field/potential
[g, voxels_in] = import_gain(path, basefile, ft_senstype(sens, 'eeg'));
if (voxels_in ~= voxels) & (nargout == 1); warning('Imported voxels from OpenMEEG process not the same as function input.'); end;
if isfield(sens,'tra'); chanmixMtx = sens.tra; else; chanmixMtx=[]; end;
lp = computeLead(g, chanmixMtx);

% Cleanup
if cleanup_flag
    rmdir(basepath,'s')
end
if (isunix & ~ismac)
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
        om_save_tri(fullfile(path, [basefile tissues{ii} '.tri']),vol.bnd(ii).pnt,vol.bnd(ii).tri);
        %savemesh([path basefile tissues{ii} '.mesh'],vol.bnd(ii).vertices/1000,vol.bnd(ii).faces-1,-vol.bnd(ii).normals);
    end
end


function Lp = computeLead(g, chanmixMtx)
if ~isempty(chanmixMtx)
    Lp = chanmixMtx*g;                 % Mchannels x (3 orientations x Nvoxels)
    Lp = Lp(diag(chanmixMtx)~=0,:);   % Take out reference sensors
else
    Lp = g;
end


function [g, voxels] = import_gain(path, basefile, eegflag)
if eegflag
    omgainfiles = dir(fullfile(path, [basefile '_eeggain*.bin']));
else
    omgainfiles = dir(fullfile(path, [basefile '_meggain*.bin']));
end
omvoxfiles = dir(fullfile(path, [basefile '_voxels*.bin']));

g=[];
voxels=[];

% join gain/voxel files
% [openmeeg calculation may have been split for memory reasons]
for ii=1:length(omgainfiles)
    g = [g om_load_full(fullfile(path, omgainfiles(ii).name),'binary')];
    voxels = [voxels;om_load_full(fullfile(path, omvoxfiles(ii).name),'binary')];
end
voxels = voxels(1:3:end,1:3);
