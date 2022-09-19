function [srf] = readsrf(filename, NoNeg)
%
%  [srf] = readsrf(filename, [NoNeg])
%
% filename: filename of .srf surface file
% NoNeg: Optional parameter for dealing with nodes that should not be
%        displayed
%        1 (default) tells the function to remove triangles containing
%                    such nodes from the mesh
%        0           leaves mesh untouched.
%
% readsrf attempts to read the BrainVoyagerQX v. 4 srf surface files
% To display:
% colormap(srf.cmap)
% trisurf(srf.triangles, srf.VX, srf.VY, srf.VZ, srf.mesh_color)
% but that would probably grind your machine to a halt, just to
% get the idea you could do
% trisurf(srf.triangles(1:100:end), srf.VX, srf.VY, srf.VZ, srf.mesh_color)
%
% BrainVoyager flags nodes that should not be displayed. The default
% behavior of this function is to remove triangles containing such nodes
% from the triangle list. If you do not like that, set NoNeg to 0.
%
% Kate Fissell 3/06
% Modified by ZSS, SSCC/NIMH/NIH
% With thanks to Hester Breman

if (nargin == 1) NoNeg = 1; end

fp = fopen(filename,'r');
if (fp == -1)
	fprintf(1,'\nError opening %s\n',filename);
	return;
end


%% read some header fields
srf.version = fread(fp,1,'float32',0,'ieee-le');
srf.reserve = fread(fp,1,'int32',0,'ieee-le');
srf.numvert = fread(fp,1,'int32',0,'ieee-le');
srf.numtri = fread(fp,1,'int32',0,'ieee-le');
srf.meshcenXYZ = fread(fp,3,'float32',0,'ieee-le');

%% print some header fields
fprintf(1,'\nVersion: %.6f',srf.version);
fprintf(1,'\nNumber of vertices: %d',srf.numvert);
fprintf(1,'\nNumber of srf.triangles: %d',srf.numtri);
fprintf(1,'\nMesh center: %.3f %.3f %.3f',srf.meshcenXYZ);


%% read vertices and normals
srf.VX = fread(fp,srf.numvert,'float32',0,'ieee-le');
srf.VY = fread(fp,srf.numvert,'float32',0,'ieee-le');
srf.VZ = fread(fp,srf.numvert,'float32',0,'ieee-le');

srf.NX = fread(fp,srf.numvert,'float32',0,'ieee-le');
srf.NY = fread(fp,srf.numvert,'float32',0,'ieee-le');
srf.NZ = fread(fp,srf.numvert,'float32',0,'ieee-le');


%% read color stuff
srf.cmap = zeros(2,3);
srf.cmap(1,:) = fread(fp,3,'float32',0,'ieee-le');
srf.alpha_convex = fread(fp,1,'float32',0,'ieee-le');
srf.cmap(2,:) = fread(fp,3,'float32',0,'ieee-le');
srf.alpha_concave = fread(fp,1,'float32',0,'ieee-le');
srf.mesh_color = fread(fp,srf.numvert,'int32',0,'ieee-le');


%% read srf.neighbors of vertices to get max number of neighbors
neighbors_fileoffset = ftell(fp);
max_neigh = 0;
big_vert = 0;
for i=1:srf.numvert
	numneigh = fread(fp,1,'int32',0,'ieee-le');
	if (numneigh > max_neigh)
		max_neigh = numneigh;
		big_vert = i;
	end
	neighbors = fread(fp,numneigh,'int32',0,'ieee-le');
end
fprintf(1,'\nVertex %d has the maximum of %d neighbors',big_vert,max_neigh);

%% rewind back to start of neighbors and read and store them
fseek(fp,neighbors_fileoffset,'bof');
srf.neighbors = zeros(srf.numvert,max_neigh+1);
for i=1:srf.numvert
	n = fread(fp,1,'int32',0,'ieee-le');
	srf.neighbors(i,1) = n;
	srf.neighbors(i,2:n+1) = fread(fp,n,'int32',0,'ieee-le');
end


%% read srf.triangles
srf.triangles = fread(fp,srf.numtri*3,'int32',0,'ieee-le');
srf.triangles = reshape(srf.triangles, [3 srf.numtri]);
srf.triangles = srf.triangles' + 1;	%% matlab uses 1 based indices so inc vertex indices

%% Remove triangles having negative colors
ipos = find (srf.mesh_color >= 0);
OK = zeros(length(srf.mesh_color),1);
OK(ipos) = 1;
N_troub = length(srf.mesh_color) - length(ipos);
if (N_troub > 0),
   fprintf(2,'Note:\nHave %d nodes that should not be rendered.\n', N_troub);
   if (NoNeg),
      fprintf(2,' Removing them from triangulation.\nSee help readsrf to change this default.\n');
      keep = ones(size(srf.triangles,1),1);
      for (i=1:1:size(srf.triangles,1)),
         if (~OK(srf.triangles(i,1)) | ~OK(srf.triangles(i,2)) | ~OK(srf.triangles(i,3)) ),
            keep(i) = 0;
         end
      end
      srf.triangles = srf.triangles(find(keep),:);
   else
      fprintf(2,' Keeping them per user''s instruction.\nSee help readsrf for options.\n');
   end
end

srf.tristrip = fread(fp,1,'int32',0,'ieee-le');
if (srf.tristrip > 0)
	srf.tristripseq = fread(fp,srf.tristrip,'int32',0,'ieee-le');
end

srf.mtcfile = fscanf(fp, '%s');
fprintf(1,'\nAssociated mtc file: ');
if (length(srf.mtcfile) > 1)
	fprintf(1,'%s',srf.mtcfile);
else
	fprintf(1,'none');
end
fprintf(1,'\n\n');


[omega cnt] = fread(fp,1,'uchar',0,'ieee-le');
if (cnt ~= 0)
	fprintf(1,'\nWarning: extra elements at end of file.');

	c = 0;
	while (cnt ~= 0)
		[omega cnt] = fread(fp,1,'uchar',0,'ieee-le');
		c = c + 1;	
        end
	fprintf(1,'\nread %d more chars.\n',c);
end

fclose(fp);
