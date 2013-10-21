function bnd = ft_import_surf (cfg)

% FT_IMPORT_SURF imports various surface tri-mesh formats into a vol 
% structure. If cfg is not specified, function will ask for parameters via 
% GUI. Currently, BrainSuite, BrainVisa and FreeSurfer meshes are 
% supported. The resulting bnd structure will be ordered from outside-in.
%
% Arguments:
% cfg.mri           = the MRI data structure, coregistered with sensor 
%                     coordinates.
% cfg.minf          = the .minf file (for BrainVisa). If not entered, the 
%                     program will attempt to find the file.
% cfg.meshes        = cell array of mesh pathnames. The order does not 
%                     matter as the meshes will automatically be reordered 
%                     from outside-in. If not specified, the user will be 
%                     prompted to provide the mesh filenames via GUI.
% cfg.nodes         = array containing number of nodes for each mesh. The 
%                     number of elements in this array should be the same 
%                     as the number of meshes. If not provided, a GUI will 
%                     request this information.
%
% Copyright (C) 2013, Daniel D.E. Wong, Sarang S. Dalal
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if ~isfield(cfg,'mri')
    error('An MRI must be provided. Please see FT_IMPORT_SURF.');
end

if ~isfield(cfg,'meshes')
    [file,path] = uigetfile({'*.dfs; *.mesh; *.surf'},'Open Surfaces','MultiSelect','on');
    if ~iscell(file); file{1} = file; end;  % In case only 1 file was selected
    for ii = 1:length(file)
        cfg.meshes{ii} = fullfile(path,file{ii});
    end
%     name = 'Order Surfaces (Outside < In)';   % Now handled automatically
%     numlines = 1;
%     defaultanswer = num2cell([1:length(file)]); for ii = 1:length(defaultanswer); defaultanswer{ii} = num2str(defaultanswer{ii}); end;
%     options.Resize='on';
%     answer = inputdlg(file,name,numlines,defaultanswer,options);
%     for ii = 1:length(answer)
%         cfg.meshes{str2num(answer{ii})} = fullfile(path,file{ii});
%         fprintf('[%d] %s\n',str2num(answer{ii}),file{ii});
%     end
end

if ~isfield(cfg,'nodes')
    lines = 1;
    prompt = cell(length(cfg.meshes),1);
    for ii = 1:length(cfg.meshes)
        [~,name,ext] = fileparts(cfg.meshes{ii});
        prompt{ii} = [name ext];
    end
    title = 'Enter Max Nodes';
    def = cell(length(cfg.meshes),1);
    for ii = 1:length(cfg.meshes)
        def{ii} = num2str(ii*1000 + 1000);
    end
    answer = inputdlg(prompt,title,lines,def);
    maxnode = zeros(1,length(cfg.meshes));
    for ii = 1:length(cfg.meshes)
        if isempty(answer{ii}); answer{ii} = def{ii}; end;
        maxnode(ii) = str2num(answer{ii});
    end
else
    maxnode = cfg.nodes;
end
assert(length(maxnode)==length(cfg.meshes),'The number of entries in cfg.nodes must be the same as the number of meshes');

% Load meshes in appropriate coordinates
for ii = 1:length(cfg.meshes)
    [~,name,ext] = fileparts(cfg.meshes{ii});
    switch(ext)
        case {'.dfs';'.DFS'}
            if ~exist('readdfs.m','file')
                error('READDFS not found. Please download from: http://brainsuite.bmap.ucla.edu/matlab_scripts/readdfs.m');
            end
            dfs = readdfs(cfg.meshes{ii});
            % Convert to MRI dimensions
            %d = sign(diag(cfg.mri.transform))';
            %dfs.vertices = dfs.vertices.*repmat(d(1:3),size(dfs.vertices,1),1);
            %for jj = 1:3; dfs.vertices(:,jj) = mod(dfs.vertices(:,jj),cfg.mri.dim(jj)); end;
            %dfs.vertices = nut_voxels2mm(dfs.vertices);
            dfs.vertices = [dfs.vertices ones(size(dfs.vertices,1),1)]*cfg.mri.transform'; dfs.vertices = dfs.vertices(:,1:3);
            bnd(ii).pnt = dfs.vertices;
            bnd(ii).tri = dfs.faces;
            bnd(ii).unit = cfg.mri.unit;
        case {'.mesh';'.MESH'}
            if ~exist('loadmesh.m','file')
                error('Please include the $BRAINVISA/brainvisa/matlab/mesh directory in the Matlab path.')
            end
            
            % Make sure transformation matrix data is loaded so we can convert units to mri
            if ~exist('minftfm','var')
                if ~isfield(cfg,'minf') % Attempt to find file
                    cfg.minf = [cfg.meshes{1} '.minf'];
                end
                
                if ~exist('minftfm','var')
                    minffid = fopen(cfg.minf);
                    if minffid == -1
                        warning('Unable to locate MINF file. This is needed to properly coregister BrainVisa meshes.');
                        minftfm = eye(4);
                    else
                        hdr=fgetl(minffid);
                        tfm_idx = strfind(hdr,'''transformations'':') + 21;
                        minftfm=sscanf(hdr(tfm_idx:end),'%f,',[4 4])';
                        fclose(minffid);
                    end
                end
            end
            
            [bnd(ii).pnt,bnd(ii).tri]=loadmesh(cfg.meshes{ii});
            bnd(ii).tri = bnd(ii).tri + 1;
            bnd(ii).unit = cfg.mri.unit;
            
            if ~isfield(cfg.mri,'transformorig')
                transformorig = eye(4);
            else
                transformorig = cfg.mri.transformorig;
            end
            
            % Convert to MRI dimensions
            bnd(ii).pnt = [bnd(ii).pnt ones(size(bnd(ii).pnt,1),1)] * minftfm' / transformorig' * cfg.mri.transform';
            bnd(ii).pnt = bnd(ii).pnt(:,1:3);
        case {'.surf','.SURF'}
            if ~exist('read_surf.m','file')
                error('Please include the $FREESURFER_HOME/matlab directory in the Matlab path.')
            end
            [bnd(ii).pnt,bnd(ii).tri]=read_surf(cfg.meshes{ii});
            bnd(ii).tri = bnd(ii).tri + 1;
            bnd(ii).unit = cfg.mri.unit;
            
            % Freesurfer mesh should be in same coord sys as original MR
            if isfield(cfg.mri,'transformorig')
                bnd(ii).pnt = [bnd(ii).pnt ones(size(bnd(ii).pnt,1),1)] / cfg.mri.transformorig' * cfg.mri.transform';
                bnd(ii).pnt = bnd(ii).pnt(:,1:3);
            end
        otherwise
            warning(['Unknown extension: ' cfg.meshes{ii} '. Skipping...']);
    end
    
    % Resample mesh
    keepratio = min(1,maxnode(ii)/size(bnd(ii).pnt,1));
    [bnd(ii).pnt, bnd(ii).tri] = meshresample(bnd(ii).pnt, bnd(ii).tri, keepratio);
    
    % Check and repair mesh
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'dup');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'isolated');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'deep');
    [bnd(ii).pnt, bnd(ii).tri] = meshcheckrepair(bnd(ii).pnt, bnd(ii).tri, 'meshfix');
end

% Order meshes from outside-in (trying to be smart here)
sort_idx = zeros(length(bnd),6);
for ii = 1:3
    max_pnt = zeros(length(bnd),1);
    min_pnt = zeros(length(bnd),1);
    for jj = 1:length(bnd)
        max_pnt(jj) = max(bnd(jj).pnt(:,ii));
        min_pnt(jj) = min(bnd(jj).pnt(:,ii));
    end
    [~,I] = sort(max_pnt,'descend');
    sort_idx(:,ii*2-1) = I;
    [~,I] = sort(min_pnt,'ascend');
    sort_idx(:,ii*2) = I;
end
bnd = bnd(floor(mode(sort_idx')));

bnd = decouplesurf(bnd);

for ii = 1:length(bnd)
    [bnd(ii).pnt, bnd(ii).tri] = surfreorient(bnd(ii).pnt, bnd(ii).tri);
    bnd(ii).tri = bnd(ii).tri(:,[3 2 1]);
end


%%
% Take care of intersecting surfaces.
%%
function bnd = decouplesurf(bnd)

for ii = 1:length(bnd)-1
    % Despite what the instructions for surfboolean says, surfaces should
    % be ordered from inside-out!!
    [newnode, newelem] = surfboolean(bnd(ii+1).pnt,bnd(ii+1).tri,'decouple',bnd(ii).pnt,bnd(ii).tri);
    bnd(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(bnd(ii+1).pnt,1);
    bnd(ii+1).pnt = newnode(newnode(:,4)==2,1:3);
end