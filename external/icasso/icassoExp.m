function sR=icassoExp(sR)
%function sR=icassoExp(sR)
%
%PURPOSE
%
%To prepare Icasso result structure for exploratory analysis, i.e.,
%to compute (dis)similarity matrix, clustering, and projection. 
%
%EXAMPLES OF BASIC USAGE
%
%First we produce an Icasso result structure...
%
%   load megdata
%   sR=icassoEst('randinit',megdata,15,'lastEig',20);
%
%The following command performs the default Icasso clustering
%procedure:
%
%   sR=icassoExp(sR);
%
%The next step would be to return results launch the visualizations
%of the results: See icassoShow, icassoViz, and icassoResult. 
%
%You can customize the Icasso procedure by using this script
%as a model but changing the optional input parameters. See also
%icassoCluster and icassoProjection. 
%
%DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%icassoExp performs two functions icassoCluster and icassoProjection:
%
%Step 1 
%
%Cluster the estimates
%
%Details:
% 1. Compute similarities (S) between the estimates
% 2. Store the similarity matrix S in field sR.cluster.similarity. 
% 3. Record the similarity function into sR.cluster.simfcn, in
%    this case, 'abscorr'. Other possibility is 'power'. You can
%    also specify an explicit similarity matrix 
% 4. Compute dissimilarities D from S using function
%    sim2dis (simply D=1-S). If you want something else, you can give another
%    function name  
% 5. Compute hierarchical clustering according to the readily computed
%    dissimilarities using average-linkage strategy. Other
%    possibilities: 'CL' (comptele link) 'SL' (single-link)
% 6. Store the clustering into field sR.cluster.partition and the
%    strategy in sR.cluster.strategy 
% 7. Compute a clustering solution validity index, R-index up to L
%    clusters (according to the partition and dissimilarities D)
%    ('rdim' sets  L equal to the (reduced) data dimension)
% 9. Store the index into field sR.cluster.index.R
% 
%   sR=icassoCluster(sR,'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
%
%Step 2
%
%Compute the coordinates for correlation graph using Curvilinear
%Component Analysis. Other possibilities: 'mmds' (principal
%coordinates), 'sammon' (Sammon's projection). See icassoProjection
%for details.
%
%Note that the dissimilarities are slightly differently scaled than in clustering:
%
%Details:
%1. Compute similarity-to-dissimilarity transformation that is
%   found to be good for the proximity preserving projection: 
%   D=sqrt(1-sR.cluster.similarity) (i.e., call function sqrtsim2dis)
%2. compute CCA on this dissimilarity matrix using 75 epochs and default
%   parameters (see icassoProjection).
%
%  sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);


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

% ver 1.21 040305

sR=icassoCluster(sR,'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);

