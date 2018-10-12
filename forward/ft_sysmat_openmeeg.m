function [dsm] = ft_sysmat_openmeeg(pos, headmodel, sens, NAflag)

% FT_SYSMAT_OPENMEEG creates a volume conduction model of the
% head using the boundary element method (BEM). This function takes
% as input the triangulated surfaces that describe the boundaries and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This function implements
%   Gramfort et al. OpenMEEG: opensource software for quasistatic
%   bioelectromagnetics. Biomedical engineering online (2010) vol. 9 (1) pp. 45
%   http://www.biomedical-engineering-online.com/content/9/1/45
%   doi:10.1186/1475-925X-9-45
% and
%   Kybic et al. Generalized head models for MEG/EEG: boundary element method
%   beyond nested volumes. Phys. Med. Biol. (2006) vol. 51 pp. 1333-1346
%   doi:10.1088/0031-9155/51/5/021
%
% This link with FieldTrip is derived from the OpenMEEG project
% with contributions from Daniel Wong and Sarang Dalal, and uses external
% command-line executables. See http://openmeeg.github.io/
%
% Use as
%   dsm = ft_sysmat_openmeeg(pos, headmodel, sens, nonadaptive_flag)
%
%
% See also FT_HEADMODEL_OPENMEEG, FT_SENSINTERP_OPENMEEG, FT_PREPARE_LEADFIELD

%$Id$

ft_hastoolbox('openmeeg', 1);  % add to path (if not yet on path)
openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries
if(~ispc) % if Linux/Mac, set number of threads
    omp_num_threads = feature('numCores');
    prefix = ['export OMP_NUM_THREADS=' num2str(omp_num_threads) ';' prefix];
end

ecog = ft_getopt(headmodel,'ecog','no');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store the current path and change folder to the temporary one
cwd = pwd;
bndom = headmodel.bnd;

workdir = fullfile(tempdir,['ft_om_' datestr(now,'ddmmyyHHMMSSFFF')]);
mkdir(workdir);

try
    % Write the triangulations to file, named after tissue type.
    % OpenMEEG v2.3 and up internally adjusts the convention for surface
    % normals, but OpenMEEG v2.2 expects surface normals to point inwards;
    % this checks and corrects if needed
    bndfile = fullfile(workdir,strcat(headmodel.tissue,'.tri'));
    for ii=1:length(bndom)
        om_save_tri(bndfile{ii}, bndom(ii).pos, bndom(ii).tri);
    end
    
    condfile  = fullfile(workdir,'om.cond');
    geomfile  = fullfile(workdir,'om.geom');
    dipfile = fullfile(workdir,'dip.txt');
    dsmfile = fullfile(workdir,'dsm.bin'); % extension added directly in om_assemble function call below
    
    % write conductivity and mesh files
    bndlabel = {};
    for i=1:length(bndom)
        [dum,bndlabel{i}] = fileparts(bndfile{i});
    end
    
    om_write_geom(geomfile,bndfile,bndlabel);
    om_write_cond(condfile,headmodel.cond,bndlabel);

    % handle dipole file
    ndip = size(pos,1);
    pos = [kron(pos,ones(3,1)),kron(ones(ndip,1),eye(3))]; % save pos with each 3D orientation
    om_save_full(pos,dipfile,'ascii');
        
    omp_num_threads = feature('numCores');
    
    if(strcmp(NAflag,'yes'))
        str = ' -DSMNA';
    else
        str = ' -DSM';
    end

    om_status = system([prefix 'om_assemble' str ' ' geomfile ' ' condfile ' ' dipfile ' ' dsmfile]);
    
    if(om_status ~= 0) % status = 0 if successful
        ft_error(['Aborting OpenMEEG pipeline due to above error.']);
    end
    
    dsm = om_load_full(dsmfile,'binary');
catch
    rethrow(lasterror)
end

try
    rmdir(workdir,'s'); % remove workdir with intermediate files
catch
    disp(lasterr);
    rmdir(workdir,'s'); % remove workdir with intermediate files
    ft_error('an error occurred while running OpenMEEG');
end
