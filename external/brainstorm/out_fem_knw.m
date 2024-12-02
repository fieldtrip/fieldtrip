function out_fem_knw(FemMat, CondTensor, OutputFile)
% OUT_FEM_KNW Write Cauchy conductivity file (.knw)
%
% INPUT: 
%    - TessMat    : Brainstorm FEM head model mesh (tetrahedral or hexahedral)
%    - CondTensor : [nLayers x 6] or [xx, yy, zz, xy, yz, zx]
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

% Indices of all elements
index_elem = (1:size(FemMat.Tissue, 1))'; 

% Anisotropic: Tensor value given for each element
if (length(CondTensor) == length(index_elem))
    values = [index_elem, CondTensor];
% Isotropic: Same value for all the elements of a give tissues, or same value for all the elements
else
    if (numel(CondTensor) == 6)
        CondTensor = repmat(CondTensor(:)', length(FemMat.TissueLabels), 1);
    else
        CondTensor = reshape(CondTensor,[],6);
    end
    % Check dimensions
    if size(CondTensor,1) ~= length(FemMat.TissueLabels)
        error(['Input mismatch: ' num2str(size(CondTensor,1)) ' conductivity tensors for  '  num2str(length(FemMat.TissueLabels)) ' tissues.']);
    end
    values = [index_elem, CondTensor(FemMat.Tissue,:)];  
end

% Open file
[fid, message] = fopen(OutputFile, 'wt');
if (fid < 0)
    error(['Could not create file: ' message]);
end
% Write file
% Note: the spaces between values are important
fprintf(fid, 'BOI - TENSORVALUEFILE\n');
fprintf(fid, '========================================================\n');
fprintf(fid, '========================================================\n');
fprintf(fid,'BOI - TENSOR\n');
fprintf(fid,'         %d   %1.5f       %1.5f       %1.5f\n             %1.5f       %1.5f       %1.5f\n', values');
fprintf(fid,'EOI - TENSOR\n');
fprintf(fid,'========================================================\n');
fprintf(fid,'========================================================\n');
fprintf(fid,'EOI - TENSORVALUEFILE\n');
% Close file
fclose(fid);


