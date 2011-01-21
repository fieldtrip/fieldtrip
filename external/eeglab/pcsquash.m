% pcsquash() - compress data using Principal Component Analysis (PCA)
%              into a principal component subspace.  To project back 
%              into the original channel space, use pcexpand()
%
% Usage: 
%        >> [eigenvectors,eigenvalues] = pcsquash(data,ncomps);
%        >> [eigenvectors,eigenvalues,compressed,datamean] ...
%                                                    = pcsquash(data,ncomps);
%
% Inputs:
%   data    = (chans,frames) each row is a channel, each column a time point
%   ncomps  = numbers of components to retain
%
% Outputs: 
%   eigenvectors = square matrix of (column) eigenvectors 
%   eigenvalues  = vector of associated eigenvalues 
%   compressed   = data compressed into space of the ncomps eigenvectors
%                  with largest eigenvalues (ncomps,frames)
%                  Note that >> compressed = eigenvectors(:,1:ncomps)'*data;
%   datamean     = input data channel (row) means (used internally)
%
% Author: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 6-97 
%
% See also: pcexpand(), svd()

% Copyright (C) 2000 Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, 
% scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 01-25-02 reformated help & license, added links -ad 

function [EigenVectors,EigenValues,Compressed,Datamean]=pcsquash(matrix,ncomps)

if nargin < 1 
   help pcsquash
   return
end
if nargin < 2
   ncomps = 0;
end
if ncomps == 0
  ncomps = size(matrix,1);
end
if ncomps < 1
   help pcsquash
   return
end

data = matrix';                    % transpose data
[n,p]=size(data);                  % now p chans,n time points
if ncomps > p
   fprintf('pcsquash(): components must be <= number of data rows (%d).\n',p);
   return;
end

Datamean = mean(data,1);  % remove column (channel) means
data = data-ones(n,1)*Datamean;    % remove column (channel) means
out=data'*data/n;
[V,D] = eig(out);                  % get eigenvectors/eigenvalues
diag(D);
[eigenval,index] = sort(diag(D));
index=rot90(rot90(index));
EigenValues=rot90(rot90(eigenval))';
EigenVectors=V(:,index);

if nargout >= 3
   Compressed = EigenVectors(:,1:ncomps)'*data';
end
