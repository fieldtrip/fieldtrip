function srf = readBESAsrf(filename)

% This method reads surface files from BESA Research/MRI. For more
% information on the srf file format also see the file format description
% provided by the BrainVoyager support site:
% http://support.brainvoyager.com
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         elp file else only the name of the file should be specified. 
%
% Return:
%     srf
%         Cell array containing info on surface points stored in file.
%          
% Copyright (C) 2013, BESA GmbH
%
% File name: readBESAsrf.m
%
% Author: Robert Spangler
% Created: 2013-11-26

% Check default parameters
% no default parameters so far...

% Ouput parameter
srf = [];

% Try to open file.
FileID = fopen(filename, 'r');
if(FileID < 0)
	printf('Error opening file.\n');
	return
end

% Read version number
srf.VerNo = fread(FileID, 1, 'float');

% Read reserved bytes, must be '0'
srf.Reserved = fread(FileID, 1, 'int');

% Read number of vertices and triangles 
srf.NoVertices = fread(FileID, 1, 'int');
srf.NoTriangles = fread(FileID, 1, 'int');

% Read mesh center (X,Y,Z), default (128.0, 128.0, 128.0)
srf.MeshCenterCoord = fread(FileID, 3, 'float');

% Read sequence of (X,Y,Z) coordinates of all vertices
srf.CoordsVertices = [];
% 1. Read NrVertices X coordinates
srf.CoordsVertices(:,1) = fread(FileID, srf.NoVertices, 'float'); 
% 2. Read NrVertices Y coordinates
srf.CoordsVertices(:,2) = fread(FileID, srf.NoVertices, 'float');
% 3. Read NrVertices Z coordinates
srf.CoordsVertices(:,3) = fread(FileID, srf.NoVertices, 'float');

% Read normal (X,Y,Z) components of all vertices
% Note: For historical reasons, normals do not point outward, but inward.
%       Multiply each normal by "-1.0" for custom rendering.    
srf.ComponentsVertices = [];
% 1. Read NrVertices X components
srf.ComponentsVertices(:,1) = fread(FileID, srf.NoVertices, 'float'); 
% 2. Read NrVertices Y components
srf.ComponentsVertices(:,2) = fread(FileID, srf.NoVertices, 'float'); 
% 3. Read NrVertices Z components
srf.ComponentsVertices(:,3) = fread(FileID, srf.NoVertices, 'float'); 

% Read convex curvature color (R,G,B,A), default (0.322, 0.733, 0.98, 1.0)
srf.ConvCurvCol = fread(FileID, 4, 'float');

% Read concave curvature color (R,G,B,A), default (0.1, 0.24, 0.32, 1.0)
srf.ConcCurvCol = fread(FileID, 4, 'float');

% Read mesh color (sequence of color indices, one for each vertex)
srf.MeshCol = fread(FileID, srf.NoVertices, 'int32');

% Read nearest neighbour data for each vertex
srf.NrNeighbours = [];
srf.Neighbours = []; 
for CurrVertex=1:srf.NoVertices
    srf.NrNeighbours(CurrVertex) = fread(FileID, 1, 'int32');
    CurrNeighs = fread(FileID, srf.NrNeighbours(CurrVertex), 'int32');
    srf.Neighbours = [srf.Neighbours; CurrNeighs];
end

% Read sequence of three indices to define vertices of each triangle
srf.Triangles = [];
srf.Triangles = fread(FileID, srf.NoTriangles * 3, 'int32');
srf.Triangles = reshape(srf.Triangles, 3, [])';

% Read number of triangle strip elements
srf.NoTriangleStripElems = fread(FileID, [1 1], 'int');

% Read sequence of strip elements
if srf.NoTriangleStripElems > 0
    srf.SequenceStripElements = fread(...
        FileID, [1 NrOfTriangleStripElements], 'int32');
end;

% Read name of MTC file
srf.MTCfile = []; % FIXME

% Close file.
fclose(FileID);
end
