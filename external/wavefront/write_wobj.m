function write_wobj(OBJ,fullfilename)
% Write objects to a Wavefront OBJ file
%
% write_wobj(OBJ,filename);
%
% OBJ struct containing:
%
% OBJ.vertices : Vertices coordinates
% OBJ.vertices_texture: Texture coordinates 
% OBJ.vertices_normal : Normal vectors
% OBJ.vertices_point  : Vertice data used for points and lines   
% OBJ.material : Parameters from external .MTL file, will contain parameters like
%           newmtl, Ka, Kd, Ks, illum, Ns, map_Ka, map_Kd, map_Ks,
%           example of an entry from the material object:
%       OBJ.material(i).type = newmtl
%       OBJ.material(i).data = 'vase_tex'
% OBJ.objects  : Cell object with all objects in the OBJ file, 
%           example of a mesh object:
%       OBJ.objects(i).type='f'               
%       OBJ.objects(i).data.vertices: [n x 3 double]
%       OBJ.objects(i).data.texture:  [n x 3 double]
%       OBJ.objects(i).data.normal:   [n x 3 double]
%
% example reading/writing,
%
%   OBJ=read_wobj('examples\example10.obj');
%   write_wobj(OBJ,'test.obj');
%
% example isosurface to obj-file,
%
%   % Load MRI scan
%   load('mri','D'); D=smooth3(squeeze(D));
%   % Make iso-surface (Mesh) of skin
%   FV=isosurface(D,1);
%   % Calculate Iso-Normals of the surface
%   N=isonormals(D,FV.vertices);
%   L=sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2)+eps;
%   N(:,1)=N(:,1)./L; N(:,2)=N(:,2)./L; N(:,3)=N(:,3)./L;
%   % Display the iso-surface
%   figure, patch(FV,'facecolor',[1 0 0],'edgecolor','none'); view(3);camlight
%   % Invert Face rotation
%   FV.faces=[FV.faces(:,3) FV.faces(:,2) FV.faces(:,1)];
%
%   % Make a material structure
%   material(1).type='newmtl';
%   material(1).data='skin';
%   material(2).type='Ka';
%   material(2).data=[0.8 0.4 0.4];
%   material(3).type='Kd';
%   material(3).data=[0.8 0.4 0.4];
%   material(4).type='Ks';
%   material(4).data=[1 1 1];
%   material(5).type='illum';
%   material(5).data=2;
%   material(6).type='Ns';
%   material(6).data=27;
%
%   % Make OBJ structure
%   clear OBJ
%   OBJ.vertices = FV.vertices;
%   OBJ.vertices_normal = N;
%   OBJ.material = material;
%   OBJ.objects(1).type='g';
%   OBJ.objects(1).data='skin';
%   OBJ.objects(2).type='usemtl';
%   OBJ.objects(2).data='skin';
%   OBJ.objects(3).type='f';
%   OBJ.objects(3).data.vertices=FV.faces;
%   OBJ.objects(3).data.normal=FV.faces;
%   write_wobj(OBJ,'skinMRI.obj');
%
% Function is written by D.Kroon University of Twente (June 2010)

if(exist('fullfilename','var')==0)
    [filename, filefolder] = uiputfile('*.obj', 'Write obj-file');
    fullfilename = [filefolder filename];
end
[filefolder,filename] = fileparts( fullfilename);

comments=cell(1,4);
comments{1}=' Produced by Matlab Write Wobj exporter ';
comments{2}='';

fid = fopen(fullfilename,'w');
write_comment(fid,comments);

if(isfield(OBJ,'material')&&~isempty(OBJ.material))
    filename_mtl=fullfile(filefolder,[filename '.mtl']);
    fprintf(fid,'mtllib %s\n',filename_mtl);
    write_MTL_file(filename_mtl,OBJ.material)
    
end

if(isfield(OBJ,'vertices')&&~isempty(OBJ.vertices))
    write_vertices(fid,OBJ.vertices,'v');
end

if(isfield(OBJ,'vertices_point')&&~isempty(OBJ.vertices_point))
    write_vertices(fid,OBJ.vertices_point,'vp');
end

if(isfield(OBJ,'vertices_normal')&&~isempty(OBJ.vertices_normal))
    write_vertices(fid,OBJ.vertices_normal,'vn');
end

if(isfield(OBJ,'vertices_texture')&&~isempty(OBJ.vertices_texture))
    write_vertices(fid,OBJ.vertices_texture,'vt');
end

for i=1:length(OBJ.objects)
    type=OBJ.objects(i).type;
    data=OBJ.objects(i).data;
    switch(type)
        case 'usemtl'
            fprintf(fid,'usemtl %s\n',data);
        case 'f'
            check1=(isfield(OBJ,'vertices_texture')&&~isempty(OBJ.vertices_texture));
            check2=(isfield(OBJ,'vertices_normal')&&~isempty(OBJ.vertices_normal));
            if(check1&&check2)
                for j=1:size(data.vertices,1)
                    fprintf(fid,'f %d/%d/%d',data.vertices(j,1),data.texture(j,1),data.normal(j,1));
                    fprintf(fid,' %d/%d/%d', data.vertices(j,2),data.texture(j,2),data.normal(j,2));
                    fprintf(fid,' %d/%d/%d\n', data.vertices(j,3),data.texture(j,3),data.normal(j,3));
                end
            elseif(check1)
                for j=1:size(data.vertices,1)
                    fprintf(fid,'f %d/%d',data.vertices(j,1),data.texture(j,1));
                    fprintf(fid,' %d/%d', data.vertices(j,2),data.texture(j,2));
                    fprintf(fid,' %d/%d\n', data.vertices(j,3),data.texture(j,3));
                end
            elseif(check2)
                for j=1:size(data.vertices,1)
                    fprintf(fid,'f %d//%d',data.vertices(j,1),data.normal(j,1));
                    fprintf(fid,' %d//%d', data.vertices(j,2),data.normal(j,2));
                    fprintf(fid,' %d//%d\n', data.vertices(j,3),data.normal(j,3));
                end
            else
                for j=1:size(data.vertices,1)
                    fprintf(fid,'f %d %d %d\n',data.vertices(j,1),data.vertices(j,2),data.vertices(j,3));
                end
            end
        otherwise
            fprintf(fid,'%s ',type);
            if(iscell(data))
                for j=1:length(data)
                    if(ischar(data{j}))
                        fprintf(fid,'%s ',data{j});
                    else
                        fprintf(fid,'%0.5g ',data{j});
                    end
                end 
            elseif(ischar(data))
                 fprintf(fid,'%s ',data);
            else
                for j=1:length(data)
                    fprintf(fid,'%0.5g ',data(j));
                end      
            end
            fprintf(fid,'\n');
    end
end
fclose(fid);

function write_MTL_file(filename,material)
fid = fopen(filename,'w');
comments=cell(1,2);
comments{1}=' Produced by Matlab Write Wobj exporter ';
comments{2}='';
write_comment(fid,comments);

for i=1:length(material)
    type=material(i).type;
    data=material(i).data;
    switch(type)
        case('newmtl')
            fprintf(fid,'%s ',type);
            fprintf(fid,'%s\n',data);
        case{'Ka','Kd','Ks'}
            fprintf(fid,'%s ',type);
            fprintf(fid,'%5.5f %5.5f %5.5f\n',data);
        case('illum')
            fprintf(fid,'%s ',type);
            fprintf(fid,'%d\n',data);
        case {'Ns','Tr','d'}
            fprintf(fid,'%s ',type);
            fprintf(fid,'%5.5f\n',data);
        otherwise
            fprintf(fid,'%s ',type);
            if(iscell(data))
                for j=1:length(data)
                    if(ischar(data{j}))
                        fprintf(fid,'%s ',data{j});
                    else
                        fprintf(fid,'%0.5g ',data{j});
                    end
                end 
            elseif(ischar(data))
                fprintf(fid,'%s ',data);
            else
                for j=1:length(data)
                    fprintf(fid,'%0.5g ',data(j));
                end      
            end
            fprintf(fid,'\n');
    end
end

comments=cell(1,2);
comments{1}='';
comments{2}=' EOF';
write_comment(fid,comments);
fclose(fid);

function write_comment(fid,comments)
for i=1:length(comments), fprintf(fid,'# %s\n',comments{i}); end

function write_vertices(fid,V,type)
switch size(V,2)
    case 1
        for i=1:size(V,1)
            fprintf(fid,'%s %5.5f\n', type, V(i,1));
        end
    case 2
        for i=1:size(V,1)
            fprintf(fid,'%s %5.5f %5.5f\n', type, V(i,1), V(i,2));
        end
    case 3
        for i=1:size(V,1)
            fprintf(fid,'%s %5.5f %5.5f %5.5f\n', type, V(i,1), V(i,2), V(i,3));
        end
    otherwise
end
switch(type)
    case 'v'
        fprintf(fid,'# %d vertices \n', size(V,1));
    case 'vt'
        fprintf(fid,'# %d texture verticies \n', size(V,1));
    case 'vn'
        fprintf(fid,'# %d normals \n', size(V,1));
    otherwise
        fprintf(fid,'# %d\n', size(V,1));
        
end




