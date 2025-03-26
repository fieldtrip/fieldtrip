function out_fem(BstFile, OutputFile, FileFormat)
% OUT_FEM Exports a Brainstorm tetrahedral/hexahedral mesh in another file format
%
% USAGE:  out_tess(BstFile, OutputFile, FileFormat=[guess])
%         out_tess(FemMat, OutputFile, FileFormat=[guess])
%
% INPUT:
%     - BstFile    : FEM head model file from the Brainstorm database
%     - FemMat     : Loaded Brainstorm FEM head model
%     - OutputFile : Full path to output filename
%     - FileFormat : 'geo', 'msh', 'dgf', ...

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

% Guess file format
if (nargin < 3) || isempty(FileFormat)
    [fPath, fBase, fExt] = bst_fileparts(OutputFile);
    switch (fExt)
        case '.msh',  FileFormat = 'MSH';
        case '.geo',  FileFormat = 'GEO';
        case '.dgf',  FileFormat = 'DGF';
        otherwise,    error(['Unsupported file extension : "' fExt '"']);
    end
end

% Load Brainstorm file
if ischar(BstFile)
    FemMat = load(BstFile);
else
    FemMat = BstFile;
    BstFile = [];
end

% Update progress bar
[fPath, fBase, fExt] = bst_fileparts(OutputFile);
bst_progress('text', 'Export surface', ['Export mesh to file "' [fBase, fExt] '"...']);

% Switch between file formats
switch upper(FileFormat)
    case 'MSH'
        out_fem_msh(FemMat, OutputFile);
    case 'GEO'
        out_fem_geo(FemMat, OutputFile);
    case 'DGF'
        out_fem_dgf(FemMat, OutputFile); 
    otherwise
        error(['Unsupported file extension : "' FileFormat '"']);
end



