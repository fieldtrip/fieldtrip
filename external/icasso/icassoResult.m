function [Iq, A, W, S, index2centrotypes]=icassoResult(sR,L)
%function [Iq, A, W, S, index]=icassoResult(sR,[L])
%
%PURPOSE 
%
%To return results of the Icasso procedure
%
%EXAMPLE OF BASIC USAGE
%
%Get stability index, source and rows of demixing matrix: 
%
%   [Iq, A, W, S]=icassoResult(sR)
%
%INPUT
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% sR  (struct) Icasso result data structure
% [L] (scalar) number of estimate-clusters, 
%      default is to set the number of estimate-clusters equal to
%      the (reduced) data dimension.
%
%OUTPUT
%
% Iq    (vector) stability index of each estimate (see details)
% A     (matrix) estimated columns of the mixing matrix (A) =
%         pinv(W) (see details) 
% W     (matrix) estimated rows of the demixing matrix (W) (see
%         details)  
% S     (matrix) estimated independent components (see details)
% index (vector) indices to the centrotypes (centroids) of each
%         estimate-cluster (see details)  
%
%DETAILS
%
%(Rows of) W correspond to centroids (actually centrotype; see help in
%function 'centrotype' for details) of the L
%estimate-clusters. These estimates represents better the "true"
%estimate than an arbitrary estimate from a single run; however,
%the resulting W does not present a strictly orthogonal base in the
%whitened space.  
%
%S is computed by W*X where X is the original data centered.
%A is computed as a pseudoinverse of W (A=pinv(W) using MATLAB notation).
%
%Iq is the heuristic quality (stability) of the estimates. Iq(2) corresponds to
%S(2,:) and W(2,:). Ideally, Iq should be 1 to each estimate. See
%function icassoStability for details. 
%
%Output argument 'index' gives indices to the centrotype estimates
%in the Icasso result data structure , e.g.,
%icassoGet(sR,'demixingmatrix',index) would return W. See icassoGet and
%icassoIdx2Centrotype. 
%
%SEE ALSO
% icassoShow icassoGet icassoStability clusterquality
%icassoIdx2Centrotype centrotype

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 johan 100105

if nargin<2 | isempty(L),
  L=icassoGet(sR,'rdim');
  disp([sprintf('\n') 'Number of estimate-clusters not given: using (reduced)' ...
	   ' data dimension.']);
end

if isempty(sR.cluster.partition)|isempty(sR.cluster.similarity),
  error('Missing clustering information.');
end

% Check if cluster level is valid
maxCluster=size(sR.cluster.partition,1);

if L<=0 | L>maxCluster,
  error('Number of requested estimate clusters out of range');
end

index2centrotypes=icassoIdx2Centrotype(sR,'partition',sR.cluster.partition(L,:));
Iq=icassoStability(sR,L,'none');
W=icassoGet(sR,'W',index2centrotypes);
A=pinv(W);

if nargout>3,
  S=icassoGet(sR,'source',index2centrotypes);  
end





	  

	  
	  