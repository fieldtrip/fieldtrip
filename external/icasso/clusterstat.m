function Stat=clusterstat(S,partition,between)
%function Stat=clusterstat(S,partition,[between])
%
%PURPOSE
%
%To compute various intra- and extra-cluster statistics.
%
%INPUT 
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% S         (matrix) NxN similarity or distance matrix between N
%             objects. Must be symmetric.  
% partition (1xN vector) partition vector (see explanation for
%             'partition vector' in function hcluster) 
% [between] (scalar) 1: compute between-cluster similarity/distance matrix, 
%                    0: (default) don't compute
% 
%OUTPUT
% 
% All fields are vectors of size 1xK, where K is the number of clusters
% 
% Stat.N(i)             the number of objects in cluster i
%
% Stat.internal.min(i)  minimum/average/max internal similarity/distance 
% Stat.internal.avg(i)  in cluster i
% Stat.internal.max(i)  
%  Note: if there is only one item in the cluster, the internal statistics
%  are set to NaN
%
% Stat.external.min(i)  minimum/average/max external similarity/distance from
% Stat.external.avg(i)  objects of cluster i to the objects of other clusters
% Stat.external.max(i)
%
%DETAILS
%
%Partition vector divides items  into clusters: partition(i) is
%the label  of cluster that item i belongs to. Cluster labels must
%be integers 1,2,...,K where K is the number of clusters.  
%
%"Internal" for cluster k refers to all pairwise
%distances/similarities between objects belonging to the same
%cluster, i.e., for objects i for which partition(i)==k. 
%
%"External" for cluster k refers to all pairwise
%distances/similarities between the members of the cluster k and
%members of the other clusters, i.e., the whole similarity matrix
%excluding the internal similarities of cluster k.
%
%If 'between' is set to 1 the function computes also the following
%KxK matrices 
%
%Stat.between.min(i,j)
%Stat.between.avg(i,j)
%Stat.between.max(i,j)
%
%These include distances/similarities  between clusters i and j
%defined as min/average/max pairwise distances between objects
%belonging to clusters i and j, respectively. 
%
%SEE ALSO
% clusterquality
% icassoStability
% icassoRindex
% rindex

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

% ver 1.2. 100105

if nargin<2,
   error('You must give at least two input arguments.');
end

if nargin<3|isempty(between)
  between=0;
end

% Number of clusters
Ncluster=max(partition);

%Initialize the struct
Stat.internal.sum(1:Ncluster,1)=NaN;
Stat.internal.min(1:Ncluster,1)=NaN;
Stat.internal.avg(1:Ncluster,1)=NaN;
Stat.internal.max(1:Ncluster,1)=NaN;
Stat.external.sum(1:Ncluster,1)=NaN;
Stat.external.min(1:Ncluster,1)=NaN;
Stat.external.avg(1:Ncluster,1)=NaN;
Stat.external.max(1:Ncluster,1)=NaN;
Stat.index = cell(1,Ncluster);

for cluster=1:Ncluster,
  thisPartition=(partition==cluster);
  Stat.index{cluster} = find(thisPartition);
  S_=S(thisPartition,thisPartition);
  Stat.N(cluster)=size(S_,1);
  S_(eye(size(S_))==1)=[];
  if ~isempty(S_),
    Stat.internal.sum(cluster)=sum(S_);
    Stat.internal.min(cluster)=min(S_);
    Stat.internal.avg(cluster)=mean(S_);
    Stat.internal.max(cluster)=max(S_);
  end 
  if Ncluster>1,
    S_=S(thisPartition,~thisPartition);
    Stat.external.sum(cluster)=sum(S_(:));
    Stat.external.min(cluster)=min(S_(:));
    Stat.external.avg(cluster)=mean(S_(:));
    Stat.external.max(cluster)=max(S_(:));
  end
end

if between,
  Stat.between.min=zeros(Ncluster,Ncluster);
  Stat.between.max=zeros(Ncluster,Ncluster);
  Stat.between.avg=zeros(Ncluster,Ncluster);

  for i=1:Ncluster,
    Pi=find(i==partition);
    for j=i+1:Ncluster,
      Pj=find(j==partition);
      d_=S(Pi,Pj); 
      Stat.between.min(i,j)=min(d_(:));
      Stat.between.avg(i,j)=mean(d_(:));		
      Stat.between.max(i,j)=max(d_(:));		;		
    end
  end  
  
  Stat.between.min=Stat.between.min+ ...
    Stat.between.min';
  Stat.between.max=Stat.between.max+ ...
      Stat.between.max';
  Stat.between.avg=Stat.between.avg+ ...
      Stat.between.avg';
end
