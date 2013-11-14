function [P,Z,order]=hcluster(D,link)
%function [P,Z,order] = hcluster(D,link)
%
%PURPOSE
%
%To perform a hierarchical agglomerative clustering on a distance
%matrix. Returns partition vectors for each level of the clustering
%hierarchy and information for visualizing the clustering hierarchy
%as a dendrogram. 
%
%INPUT 
%
% D     (MxM matrix) distance matrix  
% link  (string)     agglomeration strategy, (linkage) 
%                    'al' or 'average'  (group average link.) 
%                    'sl' or 'single',  (nearest neighbor link.)
%                    'cl' or 'complete' (furthest neighbor link.) 
%OUTPUT
% 
%  P (NxN matrix) contains the partitions on each level of the
%                 dendrogram (partition vectors) 
%  Z,order        arguments returned by som_linkage. Needed for 
%                 drawing dendrogram with function
%                 som_dendrogram. See also som_linkage.    
%
%DETAILS 
%
%A partition vector p represents division of objects into clusters,
%p(i) is the number of cluster that object i belongs to. Cluster
%numbers must be integers 1,2,...,k where k is the number of clusters. 
%
%Function hcluster clusters hierarchically N objects according to an NxN
%dissimilarity matrix D. D(i,j) is the distance between objects i
%and j. The function applies a hierarchical agglomerative
%clustering (single, complete, or group average linkage). It
%returns matrix P (of size NxN) where each row is a partition
%vector. The partition vectors present the clustering of the
%objects on each level L=1,...,N of the dendrogram. Let c=P(L,i);
%Now, c is the cluster that object i belongs to at level L. On each
%row P(L,:), cluster labels are integers 1,2,...,L.   
% 
%Note, that the cluster labels cannot be compared over
%partitions. The same cluster may appear at different level(s) but
%it does not necessarily have the same label in every partition
%vector P(L,:).  
%
%SEE ALSO
% som_linkage 
% som_dendrogram
% icassoCluster

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

% ver 1.21 johan 040305

switch lower(link),
 case 'al'
  link='average'; % Change Icasso name to SOM Toolbox name
 case 'cl'
  link='complete';% Change Icasso name to SOM Toolbox name
 case 'sl'
  link='single'; % Change Icasso name to SOM Toolbox name
end

N=size(D,1);
[Z,order]=som_linkage(ones(N,1),'linkage',link,'dist',D);
P=Z2partition(Z);

function partition=Z2partition(Z)
%
% function partition=Z2partition(Z)
%
% Recode SOM Toolbox presentation for hierachical clustering (Z)
% into partition vectors(s)

N=size(Z,1)+1;
C=zeros(N);
C(1,:)=1:N;

for i=2:N,
  C(i,:)=C(i-1,:); 
  C(i,(Z(i-1,1)==C(i,:))|(Z(i-1,2)==C(i,:)))=N-1+i;
end

for i=1:size(C,1),
  [u,tmp,newindex]=unique(C(i,:));
  C(i,:)=newindex(:)';
end

partition=C(end:-1:1,:);

