function out_fem_msh(FemMat, OutputFile)
% OUT_FEM_MSH Write tetrahedral mesh to DGF file (used by GMSH)
%
% INPUT: 
%    - TessMat    : Brainstorm FEM head model mesh (tetrahedral)
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
% Authors: Takfarinas Medani, 2015
%          Francois Tadel, 2020

% Check the format of the mesh
if (size(FemMat.Elements,2) ~= 4)
    error('This function handles only tetrahedral meshes.');
end
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
fprintf(fid, ['%s' EndLine], '$MeshFormat');
fprintf(fid, ['%s' EndLine], '2.2 0 8');
fprintf(fid, ['%s' EndLine], '$EndMeshFormat');
% Write nodes
nNodes = size(FemMat.Vertices,1);
fprintf(fid, ['%s' EndLine], '$Nodes');
fprintf(fid, ['%i' EndLine], nNodes);
fprintf(fid, ['%i %i %i %i' EndLine], [1:nNodes; FemMat.Vertices']);
fprintf(fid, ['%s' EndLine], '$EndNodes');
% Write elements
nElem = size(FemMat.Elements,1);
fprintf(fid, EndLine);
fprintf(fid, ['%s' EndLine], '$Elements');
fprintf(fid, ['%i' EndLine], nElem);
fprintf(fid, ['%i %i %i %i %i %i %i %i %i' EndLine],...
    [1:nElem; ...             % elm?number: Element index
     4 * ones(1,nElem); ...   % elm?type: 4 nodes per element (tetrahedral mesh)
     2 * ones(1,nElem); ...   % number?of?tags: Two tags per element (tissue index and unused value)
     FemMat.Tissue' - 1; ...  % Tag #1: Tissue index
     zeros(1,nElem); ...      % Tag #2: Unused
     FemMat.Elements']);      % Tetrahedral elements
fprintf(fid, ['%s' EndLine], '$EndElements');

% Close file
fclose(fid);


% % LATER: If we need to visualize data within gmsh
% % bloc physical Names
% fprintf(fid,'%s' EndLine],'$PhysicalNames');
% fprintf(fid,'%i ' EndLine],length(unique(newelem(:,6))));
% fprintf(fid,'%i %i %s' EndLine],'1 1 toto1');
% fprintf(fid,'%i %i %s' EndLine],'2 2 toto2');
% fprintf(fid,'%i %i %s' EndLine],'3 3 toto3');
% fprintf(fid,'%i %i %s' EndLine],'4 4 toto4');
% fprintf(fid,'%s' EndLine],'$EndPhysicalNames');
% 
% % $NodeData
% % http://geuz.org/gmsh/doc/texinfo/gmsh.html#SEC62
% fprintf(fid,'%s' EndLine],'$NodeData');
% fprintf(fid,'%i ' EndLine],1); % one string tag:
% fprintf(fid,'%s' EndLine],'"A scalar view"');
% fprintf(fid,'%i' EndLine],1); %one real tag:
% fprintf(fid,'%i' EndLine],0.0); %the time value (0.0)
% fprintf(fid,'%i' EndLine],3); % three integer tags:
% fprintf(fid,'%i' EndLine],0); %the time step (0; time steps always start at 0)
% fprintf(fid,'%i' EndLine],1); % 1-component (scalar) field
% fprintf(fid,'%i' EndLine],length(newnode)); % nb associated nodal values
% fprintf(fid,'%i %i ' EndLine],[nn,Vn']');
% fprintf(fid,'%s' EndLine],'$End$NodeData');





