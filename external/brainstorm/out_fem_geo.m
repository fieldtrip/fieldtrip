function out_fem_geo(FemMat, OutputFile)
% OUT_FEM_GEO Write tetrahedral mesh to Cauchy geometry .geo file (SimBio/NeuroFem)
%
% INPUT: 
%    - TessMat    : Brainstorm FEM head model mesh (tetrahedral)
%    - OutputFile : Full path to output file

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Takfarinas Medani, 2019
%          Francois Tadel, 2020

% Check the format of the mesh
if (size(FemMat.Elements,2) ~= 4)
    error('This function handles only tetrahedral meshes.');
end

% Nodes list
nNodes = size(FemMat.Vertices,1);
nodes = [FemMat.Vertices(1:2:end-1,:), FemMat.Vertices(2:2:end,:)];

% Open file
[fid, message] = fopen(OutputFile, 'w');
if (fid < 0)
    error(['Could not create file: ' message]);
end

% Write header
fprintf(fid, 'BOI - GEOMETRIEFILE\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, 'BOI - STEUERKARTE\n');
fprintf(fid, 'ANZAHL DER KNOTEN             :%d\n', nNodes);
fprintf(fid, 'ANZAHL DER ELEMENTE           :%d\n', size(FemMat.Elements,1));
fprintf(fid, 'GEOMETR. STRUKTUR - DIMENSION :      3\n');
fprintf(fid, 'EOI - STEUERKARTE\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, '===================================================================\n');
% Write nodes
fprintf(fid, 'BOI - KOORDINATENKARTE\n');
fprintf(fid, '   %1.7f %1.7f %1.7f     %1.7f %1.7f %1.7f \n', nodes');
if mod(nNodes,2)   % Write last node if odd number of nodes
    fprintf(fid, '   %1.7f %1.7f %1.7f\n', FemMat.Vertices(end,:));
end
fprintf(fid, 'EOI - KOORDINATENKARTE\n');
% Write elemets
fprintf(fid, '===================================================================\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, 'BOI - ELEMENTKNOTENKARTE\n');
fprintf(fid, '  303: %6d%6d%6d%6d\n', FemMat.Elements');
fprintf(fid, 'EOI - ELEMENTKNOTENKARTE\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, '===================================================================\n');
fprintf(fid, 'EOI - GEOMETRIEFILE\n');

% Close file
fclose(fid);

