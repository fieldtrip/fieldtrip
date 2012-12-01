function sR=icassoStruct(X)
%function sR=icassoStruct([X])
%
%PURPOSE
%
%To initiate an Icasso result data structure which is meant for
%storing and keeping organized all data, parameters and results
%when performing the Icasso procedure. 
%
%EXAMPLE OF BASIC USAGE
%
%  S=icassoStruct(X);
%
%creates an Icasso result structure in workspace variable S. Its
%fields are initially empty except field .signal that contains
%matrix X. 
%
%INPUT 
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
%[X] (dxN matrix) the original data (signals) consisting of N
%      d-dimensional vectors. Icasso centers the data (removes the 
%      sample mean from it) and stores it in field .signal. 
%      If the input argument is not given, or it is empty, 
%      the field is left empty.
%
%OUTPUT
%
% sR (struct) Icasso result structure that contains fields
%
%    .mode              (string)
%    .signal            (matrix)
%    .index             (matrix)
%    .fasticaoptions    (cell array)
%    .A                 (cell array)
%    .W                 (cell array)
%    .whiteningMatrix   (matrix)
%    .dewhiteningMatrix (matrix)
%    .cluster           (struct)
%    .projection        (struct)
%
%DETAILS
%
%The following table presents the fields of Icasso result
%structure. Icasso is a sequential procedure that is split into
%several phases (functions). The table shows the order in which 
%the fields are computed, the function that is used to change the
%parameters/results in the field, and lastly the phases that 
%the result depends on.  
%
%P=parameter that may be a explicit user input or a default parameter
%set by Icasso
%
%Phase Field                  Function          depends on field(s)
%                                                 
%(1)  .mode                   icassoEst         P   
%(1)  .signal                 icassoEst         P
%(1)  .index                  icassoEst        (ICA results)
%(1)  .fasticaoptions         icassoEst         P
%(1)  .A                      icassoEst        (ICA results)
%(1a) .W                      icassoEst        (ICA results)
%(1)  .whiteningMatrix        icassoEst        (ICA results)
%(1b) .dewhiteningMatrix      icassoEst        (ICA results)
%
%(2a) .cluster.simfcn         icassoCluster     P
%(2b) .cluster.similarity     icassoCluster     1a,1b,2a
%(2c) .cluster.s2d            icassoCluster     P
%(2d) .cluster.strategy       icassoCluster     P
%(2e) .cluster.partition      icassoCluster     2b-d
%(2f) .cluster.dendrogram     icassoCluster     2b-d
%(2g) .cluster.index.R        icassoCluster     2b,2c,2e
%
%(3a) .projection.method      icassoProjection  P
%(3b) .projection.parameters  icassoProjection  P
%(3c) .projection.coordinates icassoProjection  2b,3a-b
%
%icasso    performs the whole process with default parameters
%icassoEst performs phase 1
%icassoExp performs phases 2-3 with default parameters.
%
%(1) Data, ICA parameters, and estimation results
%
%   .mode (string) 
%     type of randomization ('bootstrap'|'randinit'|'both')
%
%   .signal (dxN matrix) 
%     the original data (signal) X (centered) where N is 
%     the number of samples and d the dimension
%
%   .index  (Mx2 matrix) 
%     the left column is the number of the estimation cycle, the
%     right one is the number of the estimate on that cycle.
%     See also function: icassoGet
%
%The following fields contain parameters and results of the ICA
%estimation using FastICA toolbox. More information can be found,
%e.g., from of function fastica in FastICA toolbox.
%
%   .fasticaoptions (cell array) 
%     contains the options that FastICA uses in estimation.
%
%   .A (cell array of matrices) 
%     contains mixing matrices from each estimation cycle
%
%   .W (cell array of matrices) 
%     contains demixing matrices from each estimation cycle
%
%   .whiteningMatrix (matrix) 
%     whitening matrix for original data (sR.signal)
%
%   .dewhiteningMatrix (matrix) 
%     dewhitening matrix for original data (sR.signal).
%
%(2) Mutual similarities and clustering
%
%Parameters and results of 
% -computing similarities S between the estimates, and
% -clustering the estimates
%are stored in field .cluster which has the following subfields:
%
%   .cluster.simfcn (string)
%     a string option for function icassoCluster
%     (icassoSimilarity); it tells how the mutual similarities
%     between estimates are computed.  
%
%   .cluster.similarity (MxM matrix) 
%     mutual similarities between estimates.
%
%   .cluster.s2d (string) 
%     before clustering and computing the clustering validity index 
%     the similarity matrix S stored in .cluster.similarity is
%     transformed into a dissimilarity matrix. This string is the
%     name of the subfunction that makes the transformation: there
%     is a function call
%        D=feval(sR.cluster.s2d,sR.cluster.similarity); 
%     inside icassoCluster. Note that the dissimilarity matrix
%     is not stored in the Icasso result data struct.  
%
%   .cluster.strategy (string) 
%     strategy that was used for hierarchical clustering which is
%     done on dissimilarities D 

%   .cluster.partition (MxM matrix) 
%     stores the partitions resulting clustering. Each row
%     partition(i,:), represents a division of M objects into K(i)
%     clusters (classes). On each row, clusters must be labeled
%     with integers 1,2,...,K(i), where K(i) is the number of
%     clusters that may be different on each row Example:
%     partition=[[1 2 3 4];[1 1 1 1];[1 1 2 2]] gives three
%     different partitions where partition(1,:) means every object
%     being in its own clusters; partition(2,:) means all objects
%     being in a single cluster, and in partition(3,:) objects 1&2
%     belong to cluster 1 and 3&4 to cluster 2.
%
%   .cluster.dendrogram.Z and .cluster.dendrogram.order
%     stores information needed for drawing dendrogram and
%     similarity matrix visualizations. More details in function
%     som_linkage 
%
%The following subfields of .cluster contain heuristic validity
%scores for the partitions in .cluster.partition. If the score is
%NaN it means that the validity has not been (or can't be)
%computed. 
%
%   .cluster.index.R (Mx1 vector) 
%     computed by subfunction rindex
%
%(3) Projection for visualization
%
%Parameters for performing the visualization projection are results
%of the projection can be found in field .projection.
%
%  .projection has the following subfields:
%
%  .projection.method (string)
%     projection method used in icassoProjection
%
%  .projection.parameters (cell array)
%     contains parameters used in icassoProjection
%
%  .coordinates (Mx2 matrix)
%     contains the coordinates of the projected estimates
%

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
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
%02111-1307, USA.

% ver 1.2 johan 100105

if nargin<1|isempty(X),
  X=[];
else
  X=remmean(X);
end

sR.mode=[]; 
sR.signal=X;
sR.index=[];
sR.fasticaoptions=[];
sR.A=[];
sR.W=[];
sR.whiteningMatrix=[];
sR.dewhiteningMatrix=[];
sR.cluster=initClusterStruct;
sR.projection=initProjectionStruct;

function cluster=initClusterStruct

cluster.simfcn=[];
cluster.similarity=[];
cluster.s2d=[];
cluster.strategy=[];
cluster.partition=[];
cluster.dendrogram.Z=[];
cluster.dendrogram.order=[];
cluster.index.R=[];

function projection=initProjectionStruct

projection.method=[];
projection.parameters=[];
projection.coordinates=[];
