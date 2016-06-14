function ri=rindex(dist,partition)
%function ri=rindex(dist,partition)
%
%PURPOSE
%
%To compute a clustering validity index called R-index to assist in
%selecting a proper number of clusters. Input argument
%'partition' defines different clustering results, partitions, of
%the same data. Clustering result that has a (local) minimum value for
%R-index among the compared ones is a "good" one.     
%
%The index is not computed (NaN) if the partition contains a
%cluster consisting of only one item.
%
%INPUT
%
% dist      (NxN matrix) matrix of pairwise dissimilarities between N objects
% partition (KxN matrix) K partition vectors of the same N objects (see
%            explanation for 'partition vector' in function hcluster) 
%
%OUTPUT
%
% ri (Kx1 vector) ri(k) is the value of R-index for clustering
%      given in partition(k,:).
%
%DETAILS
%
%R-index is a Davies-Bouldin type relative clustering validity
%index. The clustering that has the minimum value among is the best
%in terms of R-index. The index is originally defined in Levine &
%Dormay (2001),  "Resampling method for unsupervised estimation of
%cluster validity. Neural Computation, and redescribed in
%publication Himberg et al. (2004), "Validating the independent
%components of neuroimaging time-series via clustering and
%visualization". NeuroImage, 22:3(1214-1222). 13(11). 
%
%NOTE: If the partition contains clusters consisting of only one
%item, the index is especially doubtful and value NaN is assigned
%to the solution.
%
%NOTE 2: the index that is computed here is slightly
%different than in Himberg et al. (2004). In computing
%intra-cluster averages the self-similarities are not included.   
%
%SEE ALSO
% icassoRindex
% hcluster
% icassoStruct

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

% ver 1.21 030305 johan

N=size(partition,1);

% Number of clusters in each partition
Ncluster=max(partition');

disp('Computing R-index...');

for k=1:N,
  if any(hist(partition(k,:),1:Ncluster(k))==1),
    % contains one-item clusters (index very doubtful)
    ri(k,1)=NaN+1i*NaN;
  elseif Ncluster(k)==1,
    % Degenerate partition (all in the same cluster)
    ri(k,1)=NaN;
  else
    % compute cluster statistics
    s=clusterstat(dist,partition(k,:),1);
    % Compute R-index: (set diagonal to Inf to exclude self-distance!
    s.between.avg(eye(size(s.between.avg))==1)=Inf; 
    ri(k,1)=nanmean(s.internal.avg'./min(s.between.avg)); 
  end
end

