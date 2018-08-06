function saveas(this,filename,format)
% Save GIfTI object in external file format
% FORMAT saveas(this,filename,format)
% this      - GIfTI object
% filename  - name of file to be created [Default: 'untitled.vtk']
% format    - optional argument to specify encoding format, among
%             VTK (.vtk,.vtp), Collada (.dae), IDTF (.idtf). [Default: VTK]
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$


% Check filename and file format
%--------------------------------------------------------------------------
ext = '.vtk';
if nargin == 1
    filename = ['untitled' ext];
else
    if nargin == 3 && strcmpi(format,'collada')
        ext = '.dae';
    end
    if nargin == 3 && strcmpi(format,'idtf')
        ext = '.idtf';
    end
    if nargin == 3 && strncmpi(format,'vtk',3)
        format = lower(format(5:end));
        ext = '.vtk';
    end
    [p,f,e] = fileparts(filename);
    if strcmpi(e,'.gii')
        warning('Use save instead of saveas.');
        save(this,filename);
        return;
    end
    if ~ismember(lower(e),{ext})
        warning('Changing file extension from %s to %s.',e,ext);
        e = ext;
    end
    filename = fullfile(p,[f e]);
end

% Write file
%--------------------------------------------------------------------------
s = struct(this);

switch ext
    case '.dae'
        save_dae(s,filename);
    case '.idtf'
        save_idtf(s,filename);
    case {'.vtk','.vtp'}
        if nargin < 3, format = 'legacy-ascii'; end
        mvtk_write(s,filename,format);
    otherwise
        error('Unknown file format.');
end


%==========================================================================
% function save_dae(s,filename)
%==========================================================================
function save_dae(s,filename)

o = @(x) blanks(x*3);

% Split the mesh into connected components
%--------------------------------------------------------------------------
try
    C = spm_mesh_label(s.faces);
    d = [];
    for i=1:numel(unique(C))
        d(i).faces    = s.faces(C==i,:);
        u             = unique(d(i).faces);
        d(i).vertices = s.vertices(u,:);
        a             = 1:max(d(i).faces(:));
        a(u)          = 1:size(d(i).vertices,1);
        %a = sparse(1,double(u),1:1:size(d(i).vertices,1));
        d(i).faces    = a(d(i).faces);
    end
    s = d;
end

% Open file for writing
%--------------------------------------------------------------------------
fid = fopen(filename,'wt');
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% Prolog & root of the Collada XML file
%--------------------------------------------------------------------------
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<COLLADA xmlns="http://www.collada.org/2008/03/COLLADASchema" version="1.5.0">\n');

% Assets
%--------------------------------------------------------------------------
fprintf(fid,'%s<asset>\n',o(1));
fprintf(fid,'%s<contributor>\n',o(2));
fprintf(fid,'%s<author_website>%s</author_website>\n',o(3),...
    'http://www.fil.ion.ucl.ac.uk/spm/');
fprintf(fid,'%s<authoring_tool>%s</authoring_tool>\n',o(3),'SPM');
fprintf(fid,'%s</contributor>\n',o(2));
fprintf(fid,'%s<created>%s</created>\n',o(2),datestr(now,'yyyy-mm-ddTHH:MM:SSZ'));
fprintf(fid,'%s<modified>%s</modified>\n',o(2),datestr(now,'yyyy-mm-ddTHH:MM:SSZ'));
fprintf(fid,'%s<unit name="millimeter" meter="0.001"/>\n',o(2));
fprintf(fid,'%s<up_axis>Z_UP</up_axis>\n',o(2));
fprintf(fid,'%s</asset>\n',o(1));

% Image, Materials, Effects
%--------------------------------------------------------------------------
%fprintf(fid,'%s<library_images/>\n',o(1));

fprintf(fid,'%s<library_materials>\n',o(1));
for i=1:numel(s)
    fprintf(fid,'%s<material id="material%d" name="material%d">\n',o(2),i,i);
    fprintf(fid,'%s<instance_effect url="#material%d-effect"/>\n',o(3),i);
    fprintf(fid,'%s</material>\n',o(2));
end
fprintf(fid,'%s</library_materials>\n',o(1));

fprintf(fid,'%s<library_effects>\n',o(1));
for i=1:numel(s)
    fprintf(fid,'%s<effect id="material%d-effect" name="material%d-effect">\n',o(2),i,i);
    fprintf(fid,'%s<profile_COMMON>\n',o(3));
    fprintf(fid,'%s<technique sid="COMMON">\n',o(4));
    fprintf(fid,'%s<lambert>\n',o(5));
    fprintf(fid,'%s<emission>\n',o(6));
    fprintf(fid,'%s<color>%f %f %f %d</color>\n',o(7),[0 0 0 1]);
    fprintf(fid,'%s</emission>\n',o(6));
    fprintf(fid,'%s<ambient>\n',o(6));
    fprintf(fid,'%s<color>%f %f %f %d</color>\n',o(7),[0 0 0 1]);
    fprintf(fid,'%s</ambient>\n',o(6));
    fprintf(fid,'%s<diffuse>\n',o(6));
    fprintf(fid,'%s<color>%f %f %f %d</color>\n',o(7),[0.5 0.5 0.5 1]);
    fprintf(fid,'%s</diffuse>\n',o(6));
    fprintf(fid,'%s<transparent>\n',o(6));
    fprintf(fid,'%s<color>%d %d %d %d</color>\n',o(7),[1 1 1 1]);
    fprintf(fid,'%s</transparent>\n',o(6));
    fprintf(fid,'%s<transparency>\n',o(6));
    fprintf(fid,'%s<float>%f</float>\n',o(7),0);
    fprintf(fid,'%s</transparency>\n',o(6));
    fprintf(fid,'%s</lambert>\n',o(5));
    fprintf(fid,'%s</technique>\n',o(4));
    fprintf(fid,'%s</profile_COMMON>\n',o(3));
    fprintf(fid,'%s</effect>\n',o(2));
end
fprintf(fid,'%s</library_effects>\n',o(1));

% Geometry
%--------------------------------------------------------------------------
fprintf(fid,'%s<library_geometries>\n',o(1));
for i=1:numel(s)
    fprintf(fid,'%s<geometry id="shape%d" name="shape%d">\n',o(2),i,i);
    fprintf(fid,'%s<mesh>\n',o(3));
    fprintf(fid,'%s<source id="shape%d-positions">\n',o(4),i);
    fprintf(fid,'%s<float_array id="shape%d-positions-array" count="%d">',o(5),i,numel(s(i).vertices));
    fprintf(fid,'%f ',reshape(s(i).vertices',1,[]));
    fprintf(fid,'</float_array>\n');
    fprintf(fid,'%s<technique_common>\n',o(5));
    fprintf(fid,'%s<accessor count="%d" offset="0" source="#shape%d-positions-array" stride="3">\n',o(6),size(s(i).vertices,1),i);
    fprintf(fid,'%s<param name="X" type="float" />\n',o(7));
    fprintf(fid,'%s<param name="Y" type="float" />\n',o(7));
    fprintf(fid,'%s<param name="Z" type="float" />\n',o(7));
    fprintf(fid,'%s</accessor>\n',o(6));
    fprintf(fid,'%s</technique_common>\n',o(5));
    fprintf(fid,'%s</source>\n',o(4));
    fprintf(fid,'%s<vertices id="shape%d-vertices">\n',o(4),i);
    fprintf(fid,'%s<input semantic="POSITION" source="#shape%d-positions"/>\n',o(5),i);
    fprintf(fid,'%s</vertices>\n',o(4));
    fprintf(fid,'%s<triangles material="material%d" count="%d">\n',o(4),i,size(s(i).faces,1));
    fprintf(fid,'%s<input semantic="VERTEX" source="#shape%d-vertices" offset="0"/>\n',o(5),i);
    fprintf(fid,'%s<p>',o(5));
    fprintf(fid,'%d ',reshape(s(i).faces',1,[])-1);
    fprintf(fid,'</p>\n');
    fprintf(fid,'%s</triangles>\n',o(4));
    fprintf(fid,'%s</mesh>\n',o(3));
    fprintf(fid,'%s</geometry>\n',o(2));
end
fprintf(fid,'%s</library_geometries>\n',o(1));

% Scene
%--------------------------------------------------------------------------
fprintf(fid,'%s<library_visual_scenes>\n',o(1));
fprintf(fid,'%s<visual_scene id="VisualSceneNode" name="SceneNode">\n',o(2));
for i=1:numel(s)
    fprintf(fid,'%s<node id="node%d">\n',o(3),i);
    fprintf(fid,'%s<instance_geometry url="#shape%d">\n',o(4),i);
    fprintf(fid,'%s<bind_material>\n',o(5));
    fprintf(fid,'%s<technique_common>\n',o(6));
    fprintf(fid,'%s<instance_material symbol="material%d" target="#material%d"/>\n',o(7),i,i);
    fprintf(fid,'%s</technique_common>\n',o(6));
    fprintf(fid,'%s</bind_material>\n',o(5));
    fprintf(fid,'%s</instance_geometry>\n',o(4));
    fprintf(fid,'%s</node>\n',o(3));
end
fprintf(fid,'%s</visual_scene>\n',o(2));
fprintf(fid,'%s</library_visual_scenes>\n',o(1));
fprintf(fid,'%s<scene>\n',o(1));
fprintf(fid,'%s<instance_visual_scene url="#VisualSceneNode" />\n',o(2));
fprintf(fid,'%s</scene>\n',o(1));

% End of XML
%--------------------------------------------------------------------------
fprintf(fid,'</COLLADA>\n');

% Close file
%--------------------------------------------------------------------------
fclose(fid);


%==========================================================================
% function save_idtf(s,filename)
%==========================================================================
function save_idtf(s,filename)

o = @(x) blanks(x*3);

% Compute normals
%--------------------------------------------------------------------------
if ~isfield(s,'normals')
    try
        s.normals = spm_mesh_normals(...
            struct('vertices',s.vertices,'faces',s.faces),true);
    catch
        s.normals = [];
    end
end

% Split the mesh into connected components
%--------------------------------------------------------------------------
try
    C = spm_mesh_label(s.faces);
    d = [];
    try
        if size(s.cdata,2) == 1 && (any(s.cdata>1) || any(s.cdata<0))
            mi = min(s.cdata); ma = max(s.cdata);
            s.cdata = (s.cdata-mi)/ (ma-mi);
        else
        end
    end
    for i=1:numel(unique(C))
        d(i).faces    = s.faces(C==i,:);
        u             = unique(d(i).faces);
        d(i).vertices = s.vertices(u,:);
        d(i).normals  = s.normals(u,:);
        a             = 1:max(d(i).faces(:));
        a(u)          = 1:size(d(i).vertices,1);
        %a = sparse(1,double(u),1:1:size(d(i).vertices,1));
        d(i).faces    = a(d(i).faces);
        d(i).mat      = s.mat;
        try
            d(i).cdata = s.cdata(u,:);
            if size(d(i).cdata,2) == 1
                d(i).cdata = repmat(d(i).cdata,1,3);
            end
        end
    end
    s = d;
end

% Open file for writing
%--------------------------------------------------------------------------
fid = fopen(filename,'wt');
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% FILE_HEADER
%--------------------------------------------------------------------------
fprintf(fid,'FILE_FORMAT "IDTF"\n');
fprintf(fid,'FORMAT_VERSION 100\n\n');

% NODES
%--------------------------------------------------------------------------
for i=1:numel(s)
    fprintf(fid,'NODE "MODEL" {\n');
    fprintf(fid,'%sNODE_NAME "%s"\n',o(1),sprintf('Mesh%04d',i));
    fprintf(fid,'%sPARENT_LIST {\n',o(1));
    fprintf(fid,'%sPARENT_COUNT %d\n',o(2),1);
    fprintf(fid,'%sPARENT %d {\n',o(2),0);
    fprintf(fid,'%sPARENT_NAME "%s"\n',o(3),'<NULL>');
    fprintf(fid,'%sPARENT_TM {\n',o(3));
    I = s(i).mat; % eye(4);
    for j=1:size(I,2)
        fprintf(fid,'%s',o(4)); fprintf(fid,'%f ',I(:,j)'); fprintf(fid,'\n');
    end
    fprintf(fid,'%s}\n',o(3));
    fprintf(fid,'%s}\n',o(2));
    fprintf(fid,'%s}\n',o(1));
    fprintf(fid,'%sRESOURCE_NAME "%s"\n',o(1),sprintf('Mesh%04d',i));
    %fprintf(fid,'%sMODEL_VISIBILITY "BOTH"\n',o(1));
    fprintf(fid,'}\n\n');
end

% NODE_RESOURCES
%--------------------------------------------------------------------------
for i=1:numel(s)
    fprintf(fid,'RESOURCE_LIST "MODEL" {\n');
    fprintf(fid,'%sRESOURCE_COUNT %d\n',o(1),1);
    fprintf(fid,'%sRESOURCE %d {\n',o(1),0);
    fprintf(fid,'%sRESOURCE_NAME "%s"\n',o(2),sprintf('Mesh%04d',i));
    fprintf(fid,'%sMODEL_TYPE "MESH"\n',o(2));
    fprintf(fid,'%sMESH {\n',o(2));
    fprintf(fid,'%sFACE_COUNT %d\n',o(3),size(s(i).faces,1));
    fprintf(fid,'%sMODEL_POSITION_COUNT %d\n',o(3),size(s(i).vertices,1));
    fprintf(fid,'%sMODEL_NORMAL_COUNT %d\n',o(3),size(s(i).normals,1));
    if ~isfield(s(i),'cdata') || isempty(s(i).cdata)
        c = 0;
    else
        c = size(s(i).cdata,1);
    end
    fprintf(fid,'%sMODEL_DIFFUSE_COLOR_COUNT %d\n',o(3),c);
    fprintf(fid,'%sMODEL_SPECULAR_COLOR_COUNT %d\n',o(3),0);
    fprintf(fid,'%sMODEL_TEXTURE_COORD_COUNT %d\n',o(3),0);
    fprintf(fid,'%sMODEL_BONE_COUNT %d\n',o(3),0);
    fprintf(fid,'%sMODEL_SHADING_COUNT %d\n',o(3),1);
    fprintf(fid,'%sMODEL_SHADING_DESCRIPTION_LIST {\n',o(3));
    fprintf(fid,'%sSHADING_DESCRIPTION %d {\n',o(4),0);
    fprintf(fid,'%sTEXTURE_LAYER_COUNT %d\n',o(5),0);
    fprintf(fid,'%sSHADER_ID %d\n',o(5),0);
    fprintf(fid,'%s}\n',o(4));
    fprintf(fid,'%s}\n',o(3));
    
    fprintf(fid,'%sMESH_FACE_POSITION_LIST {\n',o(3));
    fprintf(fid,'%d %d %d\n',s(i).faces'-1);
    fprintf(fid,'%s}\n',o(3));
    
    fprintf(fid,'%sMESH_FACE_NORMAL_LIST {\n',o(3));
    fprintf(fid,'%d %d %d\n',s(i).faces'-1);
    fprintf(fid,'%s}\n',o(3));
    
    fprintf(fid,'%sMESH_FACE_SHADING_LIST {\n',o(3));
    fprintf(fid,'%d\n',zeros(size(s(i).faces,1),1));
    fprintf(fid,'%s}\n',o(3));
    
    if c
        fprintf(fid,'%sMESH_FACE_DIFFUSE_COLOR_LIST {\n',o(3));
        fprintf(fid,'%d %d %d\n',s(i).faces'-1);
        fprintf(fid,'%s}\n',o(3));
    end
    
    fprintf(fid,'%sMODEL_POSITION_LIST {\n',o(3));
    fprintf(fid,'%f %f %f\n',s(i).vertices');
    fprintf(fid,'%s}\n',o(3));
    
    fprintf(fid,'%sMODEL_NORMAL_LIST {\n',o(3));
    fprintf(fid,'%f %f %f\n',s(i).normals');
    fprintf(fid,'%s}\n',o(3));
    
    if c
        fprintf(fid,'%sMODEL_DIFFUSE_COLOR_LIST {\n',o(3));
        fprintf(fid,'%f %f %f\n',s(i).cdata');
        fprintf(fid,'%s}\n',o(3));
    end
        
    fprintf(fid,'%s}\n',o(2));
    fprintf(fid,'%s}\n',o(1));
    fprintf(fid,'}\n');
end

% Close file
%--------------------------------------------------------------------------
fclose(fid);
