function S = sim_struct(type)
%
% create an empty data structure
% to set up a simnibs simulation
%
% S = sim_struct(type)
% 
% A. Thielscher, 2019

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2019 Axel Thielscher, Andre Antunes, Guilherme Saturnino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.




validtypes={'SESSION', 'TMSLIST', 'TDCSLIST', 'POSITION', ...
    'ELECTRODE', 'COND', 'LIST','FIDUCIALS', 'LEADFIELD', 'TDCSLEADFIELD'};

if ~any(strcmp(type,validtypes))
    disp(validtypes)
    error('structure type has to be one of the above')
end

S.type=type;

switch S.type
       
    case 'SESSION'
       S.org=[];  % used for parsing neuronavigation data; not stored permanently; optional
       S.fname=''; % string; points towards file containing neuronavigation data; optional
       S.date=''; % string; optional
       S.poslist={}; % can be 'TMSLIST' or 'TDCSLIST'
       S.fnamehead=''; % same as ${subID}.msh created by mri2mesh or headreco
       S.subpath = ''; % path to the 'm2m_{subID}' folder created by mri2mesh or headreco (OPTIONAL, filled from fnamehead)
       S.pathfem='';   % path to save the results (OPTIONAL, filled from fnamehead)
       S.fname_tensor = ''; % file name of the diffusion tensors (OPTIONAL, filled from fnamehead)
       S.eeg_cap = ''; % file name of the CSV file with electrode positions; (OPTIONAL, filled from fnamehead)
       S.open_in_gmsh=false; % Open result in gmsh when done
       S.map_to_surf=false; % map results on individual surface (read out in middle of GM sheet)
       S.map_to_fsavg=false; % map results further onto fsaverage template
       S.map_to_vol=false; % write fields as nifti
       S.map_to_MNI=false; % write fields as nifti in MNI space
       S.tissues_in_niftis = '';  % determines for which tissues the fields will be written
                                  % to the NifTI volumes (for map_to_vol and map_to_MNI);
                                  % either 'all' or a list of tags (standard: 2 for GM)
       S.fields= 'eE'; % type of results saved in final mesh; a string containing a combination of:
                       % e (electric field strength) E (electric field vector) j (current density strength) 
                       % J (current density vector) v (electric potential) D (dA/dt vector) s (conductivity)
       S.fiducials = sim_struct('FIDUCIALS');
                
    case 'FIDUCIALS'
        S.Nz = [];
        S.Iz = [];
        S.LPA = [];
        S.RPA = [];
            
    case 'LIST'
        S.org=[];  % used for parsing neuronavigation data; not stored permanently; optional
        S.fname=''; % string; points towards neuronavigation data; optional
        S.name=''; % string; name of simulation, will be used as part of the names of the output files; optional
        S.cond=standard_cond;   % list of conductivities
        S.anisotropy_type = 'scalar'; % can be 'scalar' (use isotropic values), 'dir' (direct mapping),'mc' (mean conductivity from direct mapping),'vn' (volume normalized); optional
        S.aniso_maxratio = 10; % maximal ratio between largest eigenvalue and the two other eigenvalues of conductivity tensor
        S.aniso_maxcond = 2; % maximal directional conductivity in [S/m] (i.e. max eigenvalue of conductivity tensor)
        S.solver_options = ''; % Options to be used by the FEM solver (default is CG+AMG)
        
    case 'TMSLIST'
        S=sim_struct('LIST');
        S.type='TMSLIST';
        S.fnamecoil='';      % to chose from inside resources/coil_models
        S.pos=sim_struct('POSITION');    % list of coil positions
        
    case 'TDCSLIST'
        S=sim_struct('LIST');
        S.type='TDCSLIST';
        S.currents=[]; % currents in A (not mA!; given per stimulator channel, not electrode!)
        S.electrode=sim_struct('ELECTRODE'); % list of electrodes
        S.fnamefem=''; % name of final mesh containing results; optional

    case 'POSITION'
        S.name='';  % string; optional
        S.date='';  % string; used to store time and date at which position was safed in neuronavigation system; optional
        S.istrig=false; % boolean; used indicate whether position was automatically created by neuronavigation system after input trigger; optional
        S.matORG=[]; % position in cooridnate system used by neuronavigation system; optional
        S.orient=''; % orientation used by neuronavigation system; optional
        S.matsimnibs=[]; % 4x4 matrix defining the coil position and direction
        %[x y z c]
        %[0 0 0 1]
        % "x", "y" and "z" are column vectors defining the directions
        % of the coil. "c" is the center of the coil
        % The "y" direction is the direction the coil handle
        % The "z" direction is normal to the coil, pointing towards the head
        % The "x" direction is orthogonal to both
        % "c" is the center of the coil
        S.didt=1e6; % corresponding to 1 A/mue_s
        S.fnamefem=''; % name of final mesh containing results; optional
        %The 3 variables above offer an alternative way to define a coil
        %position. If matsimnibs is defined, they are ignored.
        S.centre = []; % center of the coil (will be projected to the head)
        S.pos_ydir = []; % Reference direction for the prolongation of the coil handle.
        S.distance = 4; % Distance from the coil to the head


    case 'ELECTRODE'
        S.name = '';
        S.definition = 'plane'; % how the electrode's position is evaluated: plane (2D) or conf (vertices - 3D); default: plane
        
        % for definition = 'plane'
        S.shape=''; % 'rect' or 'ellipse' or 'custom' for the general polygon
        S.centre = []; % center of the electrode in conf coordinate
        S.dimensions = []; % 1x2 vector with electrode size in X and Y for rectangle / ellipse shapes in [mm]
        S.pos_ydir = []; % second position used to define electrode orientation
        % The electrode's y axis is defined as (pos_ydir - centre)
        
        % for definition = 'plane' and 'conf'
        S.vertices = []; % array of vertices (nx2) for (S.definition=='plane' AND S.shape='custom')
                         % or (nx3) for S.definition=='conf'
        
        % common fields:
        S.thickness = [];  % can have up to 3 arguments. 1st argument is the lower part of the electrode
                           % 2nd arrgument is the rubber; 3rd is the upper layer
        S.channelnr = [];  % channel to which electrode is connected
        
        S.dimensions_sponge = []; % 1x2 vector with sponge size in X and Y for rectangle / ellipse shapes in [mm]; has to be larger than S.dimensions; optional
        S.holes=struct([]); % definition of a hole (reuses ELECTRODE structure, but thickness, channelnr, dimensions_sponge are not set); optional
        S.plug=struct([]);  % definition of a hole (reuses ELECTRODE structure, but thickness, channelnr, dimensions_sponge are not set); optional

    case 'COND'
        S.name = '';
        S.value = '';
        S.descrip = '';
        
    case 'LEADFIELD' % Superclass for leadfields
        S.type='LEADFIELD';
        S.date=''; % string; optional
        S.fnamehead=''; % same as ${subID}.msh created by mri2mesh or headreco
        S.subpath = ''; % path to the 'm2m_{subID}' folder created by mri2mesh or headreco (OPTIONAL, filled from fnamehead)
        S.pathfem='';   % path to save the results (OPTIONAL, filled from fnamehead)
        S.field='E';   % Field to be stored in the leadfield. Possible options are 'E'and 'J'
        S.fname_tensor = ''; % file name of the diffusion tensors (OPTIONAL, filled from fnamehead)
        S.interpolation='middle gm'; % interpolate solution to surface(s). Valid args are '' or [] (no interp), 'middle gm' (interp to middle gm layer) [default], cell array of filenames defining surfaces to interp to 
        S.interpolation_tissue = 2; % If 'interpolation' is specified, this determines the tissue compartment(s) from which to interpolate the field. Interpolation tissues must be volumes, i.e., < 1000. Specify as 1 or [1, 2, ...].
        S.tissues = 1006; % list, tissues where to store the leadfield in addition to interpolation
        S.name=''; % string; name of simulation, will be used as part of the names of the output files; optional
        S.cond=standard_cond;   % list of conductivities
        S.anisotropy_type = 'scalar'; % can be 'scalar' (use isotropic values), 'dir' (direct mapping),'mc' (mean conductivity from direct mapping),'vn' (volume normalized); optional
        S.aniso_maxratio = 10; % maximal ratio between largest eigenvalue and the two other eigenvalues of conductivity tensor
        S.aniso_maxcond = 2; % maximal directional conductivity in [S/m] (i.e. max eigenvalue of conductivity tensor)
        S.solver_options = ''; % Options to be used by the FEM solver (default is CG+AMG)

    case 'TDCSLEADFIELD'
        S=sim_struct('LEADFIELD');
        S.type='TDCSLEADFIELD';
        S.eeg_cap = ''; % file name of the CSV file with electrode positions; (OPTIONAL, filled from fnamehead)
        S.electrode = sim_struct('ELECTRODE'); % Electrode to be used. If a single electrode, repeat it. If a list, use a different electrode for each position
        S.electrode.shape = 'ellipse';
        S.electrode.dimensions = [10, 10];
        S.electrode.thickness = 4;
        

    otherwise
        error('unknown structure type')
end
