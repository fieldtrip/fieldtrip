function bd = mne_block_diag(A,n);
%
%   function bd = mne_block_diag(A,n)
%
%   Make or extract a sparse block diagonal matrix
%
%   If A is not sparse, then returns a sparse block diagonal "bd", diagonalized from the
%   elements in "A".
%   "A" is ma x na, comprising bdn=(na/"n") blocks of submatrices.
%   Each submatrix is ma x "n", and these submatrices are
%   placed down the diagonal of the matrix.
%
%   If A is already sparse, then the operation is reversed, yielding a block
%   row matrix, where each set of n columns corresponds to a block element
%   from the block diagonal.
%
%   Routine uses NO for-loops for speed considerations.



%
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
%
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% Author: John C. Mosher 1993 - 2004
%
%
% Modifications for mne Matlab toolbox
%
%   Matti Hamalainen
%   2006
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%


if(~issparse(A)),        % then make block sparse
    [ma,na] = size(A);
    bdn = na/n;             % number of submatrices

    if(bdn - fix(bdn)),
        error('Width of matrix must be even multiple of n');
    end

    tmp = reshape([1:(ma*bdn)]',ma,bdn);
    i = zeros(ma*n,bdn);
    for iblock = 1:n,
        i((iblock-1)*ma+[1:ma],:) = tmp;
    end

    i = i(:);             % row indices foreach sparse bd


    j = [1:na];
    j = j(ones(ma,1),:);
    j = j(:);             % column indices foreach sparse bd

    bd = sparse(i,j,A(:));

else                 % already is sparse, unblock it

    [mA,na] = size(A);        % matrix always has na columns
    % how many entries in the first column?
    bdn = na/n;            % number of blocks
    ma = mA/bdn;            % rows in first block

    % blocks may themselves contain zero entries.  Build indexing as above
    tmp = reshape([1:(ma*bdn)]',ma,bdn);
    i = zeros(ma*n,bdn);
    for iblock = 1:n,
        i((iblock-1)*ma+[1:ma],:) = tmp;
    end

    i = i(:);             % row indices foreach sparse bd


    j = [0:mA:(mA*(na-1))];
    j = j(ones(ma,1),:);
    j = j(:);

    i = i + j;

    bd = full(A(i));     % column vector
    bd = reshape(bd,ma,na);    % full matrix
end

return
