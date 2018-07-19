function headmodel = ft_headmodel_openmeeg(bnd, varargin)

% FT_HEADMODEL_OPENMEEG creates a volume conduction model of the
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
%   headmodel = ft_headmodel_openmeeg(bnd, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector, conductivity of each compartment
%   tissue           = tissue labels for each compartment
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD, FT_SYSMAT_OPENMEEG,
%          FT_SENSINTERP_OPENMEEG

%$Id$

ft_hastoolbox('openmeeg', 1);  % add to path (if not yet on path)
openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries
if(~ispc) % if Linux/Mac, set number of threads
    omp_num_threads = feature('numCores');
    prefix = ['export OMP_NUM_THREADS=' num2str(omp_num_threads) ' && ' prefix];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the first part is largely shared with the dipoli and bemcp implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');
tissue          = ft_getopt(varargin, 'tissue');

% copy the boundaries from the mesh into the volume conduction model
if isfield(bnd, 'bnd')
    bnd = bnd.bnd;
end

% rename pnt into pos
bnd = fixpos(bnd);

% start with an empty volume conductor
headmodel = [];

% determine the number of compartments
numboundaries = length(bnd);

% determine the desired nesting of the compartments
order = surface_nesting(bnd, 'outsidefirst');

% rearrange boundaries and conductivities
if numel(bnd)>1
    fprintf('reordering the boundaries to: ');
    fprintf('%d ', order);
    fprintf('\n');
    % update the order of the compartments
    bnd = bnd(order);
end

if isempty(conductivity)
    ft_warning('No conductivity is declared, using default values\n')
    if numboundaries == 1
        conductivity = 1;
    elseif numboundaries == 3
        % skin/skull/brain
        conductivity = [0.33 0.0042 0.33]
    elseif numboundaries == 4
        conductivity = [0.33 0.0042 1 0.33]
    else
        ft_error(['Conductivity values are required for ' num2str(numboundaries) ' shells'])
    end
    headmodel.cond = conductivity;
else
    if numel(conductivity)~=numboundaries
        ft_error('a conductivity value should be specified for each compartment');
    end
    % update the order of the compartments
    headmodel.cond = conductivity(order);
end

% assign default tissue labels if none provided 
if(isempty(tissue))
    switch(numboundaries)
        case 3
            tissue = {'scalp','skull','brain'};
        case 4
            tissue = {'scalp','skull','CSF','brain'};
        otherwise
            tissue = strcat({'domain'},num2str((1:numboundaries)'));
    end
else
    tissue = tissue(order);
end
headmodel.tissue = tissue;


headmodel.skin_surface = 1;
headmodel.source = numboundaries;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir = fullfile(tempdir,['ft_om_' datestr(now,'ddmmyyHHMMSSFFF')]);
mkdir(workdir);

try
    % Write the triangulations to file, named after tissue type.
    % OpenMEEG v2.3 and up internally adjusts the convention for surface
    % normals, but OpenMEEG v2.2 expects surface normals to point inwards;
    % this checks and corrects if needed
    bndfile = fullfile(workdir,strcat(tissue,'.tri'));
    for ii=1:length(bnd)
        ok = checknormals(bnd(ii));
        if ~ok
            % Flip faces for openmeeg convention (inwards normals)
            bndtmp = bnd(ii);
            bndtmp.tri = fliplr(bnd(ii).tri);
            
            ok = checknormals(bndtmp);
            if(ok)
                fprintf('flipping normals for OpenMEEG''s convention\n')
                bnd(ii).tri = bndtmp.tri;
            else
                % if still not ok after flip, throw a warning
                ft_warning('neither orientation of surface normals passes check... leaving as is.')
            end
            clear bndtmp
        end
        
        om_save_tri(bndfile{ii}, bnd(ii).pos, bnd(ii).tri);
    end
    
    % retain surfaces in headmodel structure (after possible normal flip)
    headmodel.bnd = bnd;
        
     
    condfile  = fullfile(workdir,'om.cond');
    geomfile  = fullfile(workdir,'om.geom');
    hmfile    = fullfile(workdir,'hm.bin');
    hminvfile = fullfile(workdir,'hminv.bin');
    
    % write conductivity and mesh files
    bndlabel = {};
    for i=1:length(bnd);
        [~,bndlabel{i}] = fileparts(bndfile{i});
    end

    om_write_geom(geomfile,bndfile,bndlabel);
    om_write_cond(condfile,headmodel.cond,bndlabel);
    
    om_status = system([prefix 'om_assemble -HM ' geomfile ' ' condfile ' ' hmfile]);
    if(om_status ~= 0) % status = 0 if successful
        ft_error('Aborting OpenMEEG pipeline due to above error.');
    end
    
    headmodel.mat = inv(om_load_sym(hmfile,'binary'));

    rmdir(workdir,'s'); % remove workdir with intermediate files
catch
    disp(lasterr);
    rmdir(workdir,'s'); % remove workdir with intermediate files
    ft_error('an error occurred while running OpenMEEG');
end

% remember the type of volume conduction model
headmodel.type = 'openmeeg';

function ok = checknormals(bnd)
points = bnd.pos;
faces = bnd.tri;

% FIXME: this method is rigorous only for star shaped surfaces
ok = 0;
% translate to the center
org = mean(points,1);
points(:,1) = points(:,1) - org(1);
points(:,2) = points(:,2) - org(2);
points(:,3) = points(:,3) - org(3);

w = sum(solid_angle(points, faces));

if w<0 && (abs(w)-4*pi)<1000*eps
    ok = 0;
    disp('Your surface normals are outwards oriented...')
elseif w>0 && (abs(w)-4*pi)<1000*eps
    ok = 1;
    %   ft_warning('your normals are inwards oriented')
else
    ft_error('your surface probably is irregular\n')
    ok = 0;
end