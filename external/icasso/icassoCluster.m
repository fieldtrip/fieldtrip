function sR=icassoCluster(sR,varargin)
%function sR=icassoCluster(sR,['identifier1',value1,'identifier2',value2,...])
%
%PURPOSE 
%
%To cluster the ICA estimates and to compute a relative clustering validity
%index (R-index)
%
%EXAMPLE OF BASIC USAGE
%
% sR=icassoCluster(sR); 
% 
%where sR is an Icasso result structure. This applies hierarchical
%clustering using group-average linkage agglomeration strategy and
%stores the results back into workspace variable sR.
%
%INPUT 
%
% sR (struct) Icasso result data structure
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)
%
% 'simfcn' (string) 'abscorr' (default) | (matrix) 
%   Indicates how to compute similarities S between estimates i,j
%   'abscorr' S(i,j) is the absolute value of the linear
%     correlation coefficient  
%   (matrix) explicitly given similarity matrix S
%     (elements should be 0...1)  
% 'strategy' (string) 'AL' (default) | 'SL' | 'CL' 
%  Sets the clustering strategy: 
%  'AL' hierarchical group average linkage
%  'SL' hierarchical single linkage (nearest neighbor)
%  'CL' hierarchical complete linkage (furthest neighbor)
% 'L' (integer) | (string) 'rdim' (default) 
%   computes a relative clustering validity index for 2...L
%   clusters. Default string 'rdim' sets L same as the (reduced)
%   data dimension  
%'s2d' (string) (default is 'sim2dis') 
%  the name of function that is used to make the transformation
%  from similarities between IC components to dissimilarities
%  D. Default 'sim2dis' makes simply D=1-S;    
%
%OUTPUT
%
% sR (struct) updated Icasso result data structure 
%
%The function updates the fields sR.cluster.* only.
%
%DETAILS
%1. The function computes similarities between the
%estimates. See Note 1 
%2. stores the similarity matrix S into field sR.cluster.similarity
%and the method into field sR.cluster.simfcn
%2. transforms the similarities into dissimilarities (distances)
%D. See Note 2. 
%3. applies the selected clustering strategy on dissimilarities. D
%The results is a partition matrix P (of size MxM). See explanation
%in function hcluster The function stores P in field
%sR.cluster.partition (Stores also other outputs of function
%hcluster into sR.cluster.dendrogram). 
%4. computes a relative clustering validity index (see function
%rindex) for dissimilarities D and partitions P from 2 to L
%clusters. If not explicitly given (or string 'rdim' is given, L
%will be the (reduced data) dimension. Stores the result in field
%sR.cluster.index.R.          
%
%NOTE 1
%
%By default, the similarities between estimates are computed as
%linear correlation coefficients (using the demixing matrix rows and
%dewhitening matrix; see function corrw). You can also explicitly
%give any MxM similarity matrix as input. The similarities
%should be between 0...1. (M is the total number of estimates.)
%
%NOTE 2
%The clustering and clustering validity indices are computed for
%dissimilarities. Therefore, the similarities must be transformed
%into dissimilarities. By default, icassoCluster uses subfunction
%sim2dis to do this; sim2dis performs simply D=1-S. If you wish to
%make the transformation otherwise, you can set a different
%function by using input identifier - valuepair
%'s2d','myfunctionname'. icassoCluster then computes
%dissimilarities by calling 
%   D=feval('myfunctionname',sR.cluster.similarity);  
%A possible function is sqrtsim2dis 
%
%SEE ALSO
% hcluster
% som_linkage
% sqrtsim2dis
% rindex
% icassoDendrogram

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

% ver 1.21 070305 johan

% Init some variables

% total number of estimates
M=icassoGet(sR,'M');

%reduced data dimension
rdim=icassoGet(sR,'rdim');

% Set default parameters
default={'simfcn','abscorr','s2d','sim2dis','strategy','AL','L','rdim'};

%% Check optional arguments and add defaults
clusterparameters=processvarargin(varargin,default);
num_of_args=length(clusterparameters);

%% check arguments
for i=1:2:num_of_args;
  switch lower(clusterparameters{i})
    
   case 'simfcn'
    simfcn=clusterparameters{i+1};
    
    % Explicit similarity matrix?
    if isnumeric(simfcn),
      if size(simfcn,1)==M & size(simfcn,2)==M,
        sR.cluster.similarity=simfcn;
        sR.cluster.simfcn='<similarities given explicitly>';
      else 
        error('Explicitly given similarity matrix has wrong size!');
      end
    else
      % should be a string
      switch lower(simfcn)
       case 'abscorr'
        % ok
        sR.cluster.simfcn=lower(simfcn);
       otherwise
        error('''simfcn'' must be string ''abscorr'' or an MxM similarity matrix');
      end
    end
   case 's2d'
    s2dfcn=lower(clusterparameters{i+1});
    if ~ischar(s2dfcn),
      error('''s2d'' must be a string (name of a function)');
    end
    sR.cluster.s2d=s2dfcn;
   case 'l'
    L=clusterparameters{i+1};
    if isnumeric(L),
      % The user has specified max number for clusters
      
      % Check L 
      if fix(L)~=L,
        error('''L'' must be an integer.');
      elseif L<2,
        error('''L'' must be at least 2.');
      elseif L>M,
        error('''L'' cannot be more than the number of estimates.');
      end
    else
      if ~strcmp(lower(L),'rdim'),
        error('''L'' expects an integer value or ''rdim''.');
      end
      % set (reduced) data dimension
      L=icassoGet(sR,'rdim');
    end
    
    if L>100,
      warning(['R-index requested for more that 100 clusters: this can' ...
               ' be heavy...']);
    end
    
   case 'strategy'
    strategy=clusterparameters{i+1};
    if ~ischar(strategy),
      error('''strategy'' must be a string');
    end
    
    % we are case insensitive
    strategy=upper(strategy);
    sR.cluster.strategy=strategy;
    
    switch sR.cluster.strategy
     case {'AL','CL','SL'}
      ; % hierarchical clustering
     otherwise
      error(['Strategy ' strategy ' not implemented.']);
    end
   otherwise
    error(['Indentifier ' clusterparameters{i} ' not recognized.']);
  end
end

%%%% Compute similarities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(sR.cluster.simfcn)
 case '<similarities given explicitly>'
  % already handled
 case 'abscorr'
  sR.cluster.similarity=abs(corrw(icassoGet(sR,'W'),icassoGet(sR,'dewhitemat')));
  %just to make sure  
  sR.cluster.similarity(sR.cluster.similarity>1)=1; 
  sR.cluster.similarity(sR.cluster.similarity<0)=0;
end

%%%%% Convert to dissimilarities using .s2d

D=feval(sR.cluster.s2d, sR.cluster.similarity);

%%%% Make partition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sR.cluster.partition,sR.cluster.dendrogram.Z,sR.cluster.dendrogram.order]=...
    hcluster(D,sR.cluster.strategy);

%%%%% Compute cluster validity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init R
sR.cluster.index.R=ones(M,1)*NaN;
% compute

sR.cluster.index.R(1:L,1)=rindex(D,sR.cluster.partition(1:L,:));  


  
