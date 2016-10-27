function OBJ=read_wobj(fullfilename)
% Read the objects from a Wavefront OBJ file
%
% OBJ=read_wobj(filename);
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
% Example,
%   OBJ=read_wobj('examples\example10.obj');
%   FV.vertices=OBJ.vertices;
%   FV.faces=OBJ.objects(3).data.vertices;
%   figure, patch(FV,'facecolor',[1 0 0]); camlight
%
% Function is written by D.Kroon University of Twente (June 2010)

verbose=true;

if(exist('fullfilename','var')==0)
    [filename, filefolder] = uigetfile('*.obj', 'Read obj-file');
    fullfilename = [filefolder filename];
end
filefolder = fileparts( fullfilename);
if(verbose),disp(['Reading Object file : ' fullfilename]); end


% Read the DI3D OBJ textfile to a cell array
file_words = file2cellarray( fullfilename);
% Remove empty cells, merge lines split by "\" and convert strings with values to double
[ftype fdata]= fixlines(file_words);

% Vertex data
vertices=[]; nv=0;
vertices_texture=[]; nvt=0;
vertices_point=[]; nvp=0;
vertices_normal=[]; nvn=0;
material=[];

% Surface data
no=0;

% Loop through the Wavefront object file
for iline=1:length(ftype)
    if(mod(iline,10000)==0),
        if(verbose),disp(['Lines processed : ' num2str(iline)]); end
    end
    
    type=ftype{iline}; data=fdata{iline};
    
    % Switch on data type line
    switch(type)
        case{'mtllib'}
            if(iscell(data))
                datanew=[];
                for i=1:length(data)
                    datanew=[datanew data{i}];
                    if(i<length(data)), datanew=[datanew ' ']; end
                end
                data=datanew;
            end
            
            filename_mtl=fullfile(filefolder,data);
            material=readmtl(filename_mtl,verbose);
        case('v') % vertices
            nv=nv+1;
            if(length(data)==3)
                % Reserve block of memory
                if(mod(nv,10000)==1), vertices(nv+1:nv+10001,1:3)=0; end
                % Add to vertices list X Y Z
                vertices(nv,1:3)=data;
            else
                % Reserve block of memory
                if(mod(nv,10000)==1), vertices(nv+1:nv+10001,1:4)=0; end
                % Add to vertices list X Y Z W
                vertices(nv,1:4)=data;
            end
        case('vp')
            % Specifies a point in the parameter space of curve or surface
            nvp=nvp+1;
            if(length(data)==1)
                % Reserve block of memory
                if(mod(nvp,10000)==1), vertices_point(nvp+1:nvp+10001,1)=0; end
                % Add to vertices point list U
                vertices_point(nvp,1)=data;
            elseif(length(data)==2)
                % Reserve block of memory
                if(mod(nvp,10000)==1), vertices_point(nvp+1:nvp+10001,1:2)=0; end
                % Add to vertices point list U V
                vertices_point(nvp,1:2)=data;
            else
                % Reserve block of memory
                if(mod(nvp,10000)==1), vertices_point(nvp+1:nvp+10001,1:3)=0; end
                % Add to vertices point list U V W
                vertices_point(nvp,1:3)=data;
            end
        case('vn')
            % A normal vector
            nvn=nvn+1; if(mod(nvn,10000)==1),  vertices_normal(nvn+1:nvn+10001,1:3)=0; end
            % Add to vertices list I J K
            vertices_normal(nvn,1:3)=data;
        case('vt')
            % Vertices Texture Coordinate in photo
            % U V W
            nvt=nvt+1;
            if(length(data)==1)
                % Reserve block of memory
                if(mod(nvt,10000)==1), vertices_texture(nvt+1:nvt+10001,1)=0; end
                % Add to vertices texture list U
                vertices_texture(nvt,1)=data;
            elseif(length(data)==2)
                % Reserve block of memory
                if(mod(nvt,10000)==1), vertices_texture(nvt+1:nvt+10001,1:2)=0; end
                % Add to vertices texture list U V
                vertices_texture(nvt,1:2)=data;
            else
                % Reserve block of memory
                if(mod(nvt,10000)==1), vertices_texture(nvt+1:nvt+10001,1:3)=0; end
                % Add to vertices texture list U V W
                vertices_texture(nvt,1:3)=data;
            end
        case('l')
            no=no+1; if(mod(no,10000)==1), objects(no+10001).data=0; end
            array_vertices=[];
            array_texture=[];
            for i=1:length(data),
                switch class(data)
                    case 'cell'
                        tvals=str2double(stringsplit(data{i},'/'));
                    case 'string'
                        tvals=str2double(stringsplit(data,'/'));
                    otherwise
                        tvals=data(i);
                end
                val=tvals(1);
                if(val<0), val=val+1+nv; end
                array_vertices(i)=val;
                if(length(tvals)>1),
                    val=tvals(2);
                    if(val<0), val=val+1+nvt; end
                    array_texture(i)=val;
                end
            end
            objects(no).type='l';
            objects(no).data.vertices=array_vertices;
            objects(no).data.texture=array_texture;
        case('f')
            no=no+1; if(mod(no,10000)==1), objects(no+10001).data=0; end
            array_vertices=[];
            array_texture=[];
            array_normal=[];
            for i=1:length(data);
                switch class(data)
                    case 'cell'
                        tvals=str2double(stringsplit(data{i},'/'));
                    case 'string'
                        tvals=str2double(stringsplit(data,'/'));
                    otherwise
                        tvals=data(i);
                end
                val=tvals(1);
                
                if(val<0), val=val+1+nv; end
                array_vertices(i)=val;
                if(length(tvals)>1),
                    if(isfinite(tvals(2)))
                        val=tvals(2);
                        if(val<0), val=val+1+nvt; end
                        array_texture(i)=val;
                    end
                end
                if(length(tvals)>2),
                    val=tvals(3);
                    if(val<0), val=val+1+nvn; end
                    array_normal(i)=val;
                end
            end
            
            % A face of more than 3 indices is always split into
            % multiple faces of only 3 indices.
            objects(no).type='f';
            findex=1:min (3,length(array_vertices));
            objects(no).data = [];
            objects(no).data.vertices=array_vertices(findex);
            if(~isempty(array_texture)),objects(no).data.texture=array_texture(findex); end
            if(~isempty(array_normal)),objects(no).data.normal=array_normal(findex); end
            for i=1:length(array_vertices)-3;
                no=no+1; if(mod(no,10000)==1), objects(no+10001).data=0; end
                findex=[1 2+i 3+i];
                findex(findex>length(array_vertices))=findex(findex>length(array_vertices))-length(array_vertices);
                objects(no).type='f';
                objects(no).data.vertices=array_vertices(findex);
                if(~isempty(array_texture)),objects(no).data.texture=array_texture(findex); end
                if(~isempty(array_normal)),objects(no).data.normal=array_normal(findex); end
            end
        case{'#','$'}
            % Comment
            tline='  %'; 
            if(iscell(data))
                for i=1:length(data), tline=[tline ' ' data{i}]; end
            else
                tline=[tline data];
            end
            if(verbose), disp(tline); end
        case{''}
        otherwise
            no=no+1;
            if(mod(no,10000)==1), objects(no+10001).data=0; end
            objects(no).type=type;
            objects(no).data=data;
    end
end

% Initialize new object list, which will contain the "collapsed" objects
objects2(no).data=0;

index=0;

i=0;
while (i<no), i=i+1;
    type=objects(i).type;
    % First face found
    if((length(type)==1)&&(type(1)=='f'))
        % Get number of faces
        for j=i:no
            type=objects(j).type;
            if((length(type)~=1)||(type(1)~='f'))
                j=j-1; break;
            end
        end
        numfaces=(j-i)+1;
        
        index=index+1;
        objects2(index).type='f';
        % Process last face first to allocate memory
        objects2(index).data.vertices(numfaces,:)= objects(i).data.vertices;
        if(isfield(objects(i).data,'texture'))
            objects2(index).data.texture(numfaces,:) = objects(i).data.texture;
        else
            objects2(index).data.texture=[];
        end
        if(isfield(objects(i).data,'normal'))
            objects2(index).data.normal(numfaces,:)  = objects(i).data.normal;
        else
            objects2(index).data.normal=[];
        end
        % All faces to arrays
        for k=1:numfaces
            objects2(index).data.vertices(k,:)= objects(i+k-1).data.vertices;
            if(isfield(objects(i).data,'texture'))
                objects2(index).data.texture(k,:) = objects(i+k-1).data.texture;
            end
            if(isfield(objects(i).data,'normal'))
                objects2(index).data.normal(k,:)  = objects(i+k-1).data.normal;
            end
        end
        i=j;
    else
        index=index+1;
        objects2(index).type=objects(i).type;
        objects2(index).data=objects(i).data;
    end
end

% Add all data to output struct
OBJ.objects=objects2(1:index);
OBJ.material=material;
OBJ.vertices=vertices(1:nv,:);
OBJ.vertices_point=vertices_point(1:nvp,:);
OBJ.vertices_normal=vertices_normal(1:nvn,:);
OBJ.vertices_texture=vertices_texture(1:nvt,:);
if(verbose),disp('Finished Reading Object file'); end


function twords=stringsplit(tline,tchar)
% Get start and end position of all "words" separated by a char
i=find(tline(2:end-1)==tchar)+1; i_start=[1 i+1]; i_end=[i-1 length(tline)];
% Create a cell array of the words
twords=cell(1,length(i_start)); for j=1:length(i_start), twords{j}=tline(i_start(j):i_end(j)); end

function file_words=file2cellarray(filename)
% Open a DI3D OBJ textfile
fid=fopen(filename,'r');
file_text=fread(fid, inf, 'uint8=>char')';
fclose(fid);
file_lines = regexp(file_text, '\n+', 'split');
file_words = regexp(file_lines, '\s+', 'split');

function [ftype fdata]=fixlines(file_words)
ftype=cell(size(file_words));
fdata=cell(size(file_words));

iline=0; jline=0;
while(iline<length(file_words))
    iline=iline+1;
    twords=removeemptycells(file_words{iline});
    if(~isempty(twords))
        % Add next line to current line when line end with '\'
        while(strcmp(twords{end},'\')&&iline<length(file_words))
            iline=iline+1;
            twords(end)=[];
            twords=[twords removeemptycells(file_words{iline})];
        end
        % Values to double
        
        type=twords{1};
        stringdold=true;
        j=0;
        switch(type)
            case{'#','$'}
                for i=2:length(twords)
                    j=j+1; twords{j}=twords{i};                    
                end    
            otherwise    
                for i=2:length(twords)
                    str=twords{i};
                    if strcmpi(str, 'nan')
                      val=NaN;
                      stringd=false;
                    else
                      val=str2double(str);
                      stringd=~isfinite(val);
                    end
                    if(stringd)
                        j=j+1; twords{j}=str;
                    else
                        if(stringdold)
                            j=j+1; twords{j}=val;
                        else
                            twords{j}=[twords{j} val];    
                        end
                    end
                    stringdold=stringd;
                end
        end
        twords(j+1:end)=[];
        jline=jline+1;
        ftype{jline}=type;
        if(length(twords)==1), twords=twords{1}; end
        fdata{jline}=twords;
    end
end
ftype(jline+1:end)=[];
fdata(jline+1:end)=[];

function b=removeemptycells(a)
j=0; b={};
for i=1:length(a);
    if(~isempty(a{i})),j=j+1; b{j}=a{i}; end;
end

function  objects=readmtl(filename_mtl,verbose)
if(verbose),disp(['Reading Material file : ' filename_mtl]); end
file_words=file2cellarray(filename_mtl);
% Remove empty cells, merge lines split by "\" and convert strings with values to double
[ftype fdata]= fixlines(file_words);

% Surface data
objects.type(length(ftype))=0; 
objects.data(length(ftype))=0; 
no=0;
% Loop through the Wavefront object file
for iline=1:length(ftype)
  type=ftype{iline}; data=fdata{iline};
    
    % Switch on data type line
    switch(type)
        case{'#','$'}
            % Comment
            tline='  %'; 
            if(iscell(data))
                for i=1:length(data), tline=[tline ' ' data{i}]; end
            else
                tline=[tline data];
            end
            if(verbose), disp(tline); end
        case{''}
        otherwise
            no=no+1;
            if(mod(no,10000)==1), objects(no+10001).data=0; end
            objects(no).type=type;
            objects(no).data=data;
    end
end
objects=objects(1:no);
if(verbose),disp('Finished Reading Material file'); end

