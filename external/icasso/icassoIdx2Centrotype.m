function index2centrotype=icassoIdx2Centrotype(sR,mode,P)
%function index2centrotype=icassoIdx2Centrotype(sR,mode,P)
%
%PURPOSE
%
%To compute the index to centrotype of
%estimate-clusters. Centrotype is an estimate close to the centroid
%of the estimate-cluster. (see help in function centrotype for an
%exact definition). This estimate represents better the "true"
%estimate than an arbitrary estimate from a single run.         
%
%EXAMPLE OF BASIC USAGE
%
%For this example variable sR must contain complete results of 
%Icasso randomization & clustering (see functions icassoEst and
%icassoExp):
%
%   U=sR.cluster.partition(10,:); 
%   i=icassoIdx2Centrotype(sR,'partition',U); 
%   signalplot(icassoGet(sR,i));
%
%plots 10 estimated sources (ICs) each corresponding to a
%centrotype of the 10 estimate-clusters at level 10 in the dendrogram.
%
%INPUT
%
% sR   (struct) Icasso result data structure
% mode (string) 'index' | 'partition' 
% P    (vector) 1xk vector that together with mode 
%        defines the estimates of which the centrotype(s) is/are 
%        calculated
%
%        for 'partition' P is a partition vector (see
%        explanation for 'partition vector' in function hcluster)
%        Element index2centrotype(i) in output refers now to the
%        centrotype of cluster i implied by the partition vector
%        given in P. 
%
%        for 'index' P gives direct indices into the
%        estimates (see icassoGet). The centrotype is computed
%        among these estimates.
% 
%OUTPUT
% 
% index2centrotype (scalar or vector) index/indices that can be
% used, e.g., when using function icassoGet. See icassoEst for a
% note about the indexing convention.
% 
%SEE ALSO
% centrotype
% icassoGet

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

% ver 1.2 100105 johan

if isempty(sR.cluster.similarity),
  error('Similarity matrix is not computed yet!');
end

switch lower(mode)
 case 'index'
  i=centrotype(sR.cluster.similarity(P,P));
  index2centrotype=P(i);
 case 'partition'
  Ncluster=max(P);
  for cluster=1:Ncluster,
    index=find(P==cluster);
    index2centrotype(cluster,1)=icassoIdx2Centrotype(sR,'index',index);
  end
 otherwise
  error('Unknown operation mode.')
end


  