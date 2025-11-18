function out_fem_dgf(FemMat, OutputFile)
% OUT_FEM_DGF Write hexahedral mesh to MSH file (used by GMSH)
%
% INPUT: 
%    - TessMat    : Brainstorm FEM head model mesh (hexahedral)
%    - OutputFile : Full path to output file
%
% DOCUMENTATION:
%    - http://www.ensta-paristech.fr/~kielbasi/docs/gmsh.pdf

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
% Authors: Takfarinas Medani, 2019-2020
%          Francois Tadel, 2020

% Check the format of the mesh
if (size(FemMat.Elements,2) ~= 8)
    error('This function handles only hexahedral meshes.');
end
% TODO : Check if the index of the nodes and the element -1 (index on cpp)
node = FemMat.Vertices;
elem = [FemMat.Elements, FemMat.Tissue] - 1;

% Ensure end of lines are always represented by \r\n
if ispc
    EndLine = '\n';
else
    EndLine = '\r\n';
end
% Open file
[fid, message] = fopen(OutputFile, 'wt');
if (fid < 0)
    error(['Could not create file: ' message]);
end
% Write header 
fprintf(fid, ['%s' EndLine], 'DGF');
% Write nodes
fprintf(fid, ['%s' EndLine], 'Vertex');
fprintf(fid, ['%i %i %i' EndLine], node');
fprintf(fid, ['%s' EndLine], '#');
% Write elements
fprintf(fid, ['%s' EndLine], 'Cube');
fprintf(fid, ['%s' EndLine], 'parameters 1');
fprintf(fid, ['%i %i %i %i %i %i %i %i %i' EndLine], elem');
fprintf(fid, ['%s' EndLine], '#');
% Close file
fclose(fid);


