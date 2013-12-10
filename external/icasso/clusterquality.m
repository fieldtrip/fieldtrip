function [score,in,out]=clusterquality(f,S,partition)
%function [score,in,out]=clusterquality(f,S,partition)
%
%PURPOSE
%
%To compute a simple quality (compactness) index for the given
%clusters. The more dense and isolated a cluster is the bigger
%value the index gets; however, if there is only one item in a
%cluster the index is not computed (a NaN is returned for that
%cluster)    
%
%INPUT
%
% f          (string)     cluster score function: 'mean' or 'minmax'
% S          (NxN matrix) similarity matrix
% partition  (Nx1 vector) partition vector (see explanation for
%              'partition vector' in function hcluster) 
%
%OUTPUT
%
%score       (Kx1 vector) score(k) contains the index value for
%             cluster k: score(k)=in(k)-out(k) 
%in          (Kx1 vector)  see above 
%out         (Kx1 vector)  see above
%
%DETAILS
%
%For f='mean' the average extra-cluster similarity is subtracted
%from the mean intra-cluster similarity.
%
%   Iq=avg(intra-cluster similarity) - avg(extra-cluster similarity)
%
%The index is described in publication Himberg et al. (2004),
%"Validating the independent components of neuroimaging time-series
%via clustering and visualization". NeuroImage, 22:3(1214-1222). It
%is used by default in Icasso functions.
%
%(The intra-cluster similarities for cluster C mean the mutual
%similarities between the estimates in C, and the extra-cluster
%similarities for C are the similarities between the estimates in C
%and the estimates not in C.) 
%
%Note: if there is only one item in a cluster, index
%value NaN is returned for that cluster.
%
%You can also specify a more conservative but less robust criteria
%('minmax') where the maximum extra-cluster similarity is
%subtracted from the minimum intra-cluster similarity: if
%f='minmax' the maximum extra-cluster similarity is subtracted from
%the minimum intra-cluster similarity. 
%
%SEE ALSO
% clusterstat
% icassoStability

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

% ver 1.2 100105

Ncluster=max(partition);
s=clusterstat(S,partition);

% Compute score
switch lower(f)
 case 'mean'
  in=s.internal.avg;
  out=s.external.avg;
 case 'minmax'
  in=s.internal.min;
  out=s.external.max;
 otherwise
  error('Unrecognized score function.');
end

score=in-out;

