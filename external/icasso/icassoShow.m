function [Iq, A, W, S, index]=icassoShow(sR,varargin)
%function icassoShow(sR,['identifier1',value1,'indentifier2',value2,...])
%
%PURPOSE
%
%To generate explorative visualizations for Icasso
%
%EXAMPLES OF BASIC USAGE
%
%  [Iq, A, W, S]=icassoShow(sR); 
%
%shows results for as many estimate-clusters as there are (reduced)
%data dimensions. Also, return the estimates of independent
%components estimates (A,W,S) that correspond to the centroid of
%each estimate-cluster. The first output Iq contain the quality of
%the  estimates. You can rank the estimates according to this index.
%
%  icassoShow(sR,'colorlimit',[0.7 0.9],'estimate','demixing','L',9);
%
%changes the color scale
%0...0.7 (not shown), 0.7...0.9 (light red), 0.9...1 (bright red)
%and suppresses the graph lines for similarities under value 0.7 in
%general, and inside clusters that are dense (0.9...1). Shows rows
%of demixing matrix instead of sources in the estimate
%window. Aggregate results in 9 estimate-clusters.          
%  
%INPUTS
%
%sR (struct) Icasso result data structure
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)  
%
% 'L' (string) 'rdim' (default) | (integer) 
%   sets the number of estimate-clusters 'rdim' sets it equal to the
%   (reduced) data dimension 
% 'estimate' (string) 'source' (default) | 'demixing' | 'mixing' | 'off' 
%   whether to show the estimates of the
%   - independent components (sources), 
%   - rows of the demixing matrix (W), or
%   - columns of the mixing matrix (A)
%    that are associated to the centrotype of each estimate-cluster.
%   Argument 'off' suppresses the window.
% 'colorlimit' (vector) default [0.5 0.75 0.9] 
%   sets the thresholds for color of graph lines and clusters 
%   in the 2D plot; if the cluster density (the average
%   intra-cluster similarity) exceeds the highest value, the
%   cluster/lines will be bright red, if it is below the minimum,
%   the cluster is white/lines are suppressed. The rest is colored
%   with shades of red.    
% 'line' (string) 'on' (default) | 'off' 
%   whether to show the similarity graph lines in the 2D plot or not
% 'hull' (string) 'on' (default) | 'off'
%   whether to show the "cluster hulls" in the 2D plot or not
% 'graphlimit' (scalar) in 0...1 | (string) 'auto' (default) 
%   Controls the 2D plot: See function icassoGraph
% 'dense' (scalar) in 0...1 | (string) 'auto' (default) 
%   Controls the 2D plot: See function icassoGraph
% 'quality' (string) 'simple' (default) | 'detailed'
%   'simple' shows the cluster quality index, 'detailed' shows also
%   more detailed info in a separate window (figure 6)
%
%OUTPUT
%
% Iq    (vector) stability index of each estimate (see function
%         icassoResult) 
% A     (matrix) estimated columns of the mixing matrix (A) =
%         pinv(W) (see function icassoResult)
% W     (matrix) estimated rows of the demixing matrix (W) (see
%         function icassoResult) 
% S     (matrix) estimated independent components (see function
%         icassoResult) 
% index (vector) indices to the centrotypes (centroids) of each
%         estimate-cluster (see function icassoResult)  
%
%DETAILS
%
%Detailed explanation of the resulting figures can be found in help
%texts of the functions mentioned below:
%
%[Figures in brackets can be suppressed or are optional by default]
%
%Figure 1: Relative clustering quality index for different number
%of clusters and additional information: number of clusters used in
%the rest of figures (L), number of ICs, and (reduced) data
%dimension. See function icassoRindex. 
%
%Figure 2: Stability (reliability) indices  of the selected L
%estimate-clusters. See function icassoStability. 
%
%Figure 3: Correlation structure as a matrix and a dendrogram
%representation for L clusters. See function icassoDendrogram.   
%
%Figure 4: Graph of the correlations between all the estimates and
%estimate-clusters as convex hulls for L clusters
%(see function icassoGraph) 
%
%[Figure 5: IC estimates, rows of W, or columns of A (depends on user
%selection) that correspond to the centroid (actually centrotype)
%of the selected estimate-clusters. See function signalplot.] 
%
%[Figure 6: more detailed statistics on the estimate
%clusters: components of stability index Iq. See function
%icassoStability.] 
%
%SEE ALSO
% icassoGet
% icassoResult
% signalplot
% icassoStability
% icassoViz
% icasso

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

% ver 1.21 040305 johan 

if nargin<1|isempty(sR),
  error('At least one input argument expected');
end

if isempty(sR.projection.coordinates) | isempty(sR.cluster.partition) | ...
      isempty(sR.cluster.index) | isempty(sR.cluster.similarity), 
  error('Missing similarity/projection/cluster information.');
end

% initiate output args.
index2centrotypes=[]; clusterquality=[]; partition=[];

%% Set defaults and process optional input
default={'line','on','estimate','source', 'quality','simple','L',icassoGet(sR,'rdim'),...
	 'graphlimit','auto','colorlimit',[0.5 0.75 0.9],...
	 'dense','auto','hull','on'};

% initiate arguments to icassoGraph
graphArgs=[];

varargin=processvarargin(varargin,default);
num_of_args=length(varargin);

for i=1:2:num_of_args,
  id=varargin{i}; value=varargin{i+1};
  switch lower(id)
   % Check first icassoShow
   case 'quality'
    switch lower(value)
     case 'simple'
      detailedrankplot=0;
     case 'detailed'
      detailedrankplot=1;
     otherwise
      error('Option ''quality'' must be ''simple'' or ''detailed''.');
    end
   case 'estimate'
    switch lower(value)
     case {'demixing','mixing','off','source'}
      est=lower(value);
     otherwise
      error(['Option ''estimate'' must be ''source'',''demixing'',' ...
	     ' ''mixing'', or ''off''.']);
    end
   case 'l'
    if isnumeric(value);
      level=value;
    else
      switch lower(value)
       case 'rdim'
	level=icassoGet(sR,'rdim');
       otherwise
	error(['Option ''L'' must be an integer or string ''rdim''']);
      end
    end
    % submit also to icassoGraph
    graphArgs{1,end+1}=id;
    graphArgs{1,end+1}=value;
    % submit the following to icassoGraph
   case {'line','graphlimit','colorlimit','dense','hull'}
    graphArgs{1,end+1}=id;
    graphArgs{1,end+1}=value;
   otherwise
    error(['Option ''' lower(id) ''' not available.' sprintf('\n') ...
	   'Available: ''L'', ''line'',''graphlimit'',''colorlimit'',''hull'',' ...
	   '''estimate'', and ''quality''.']);
  end
end

% Check if cluster level is valid
maxCluster=size(sR.cluster.partition,1);
if level<=0 | level>maxCluster,
  error('Cluster level out of range or not specified.');
end

%%%% Get main results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Iq, A, W, S, index2centrotypes]=icassoResult(sR,level);

%%%%% Compute some cluster statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the partition

partition=sR.cluster.partition(level,:);
Ncluster=max(partition);

% cluster statistics
c=sR.cluster.similarity;
s=clusterstat(c,partition);

%%%%%%%%%%%%%%%%% Clustering validity index %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf reset;

icassoRindex(sR,level);

%%%%%%%%%%%%%%%% Ranking clusters & number of estimates %%%%%%%%%%%%%%%

figure(2); clf;

% compute & plot quality index

subplot(1,2,1);
Iq=icassoStability(sR,level,'plotindex');

% compute the rank order of estimates
clusterlabels=1:Ncluster;
[tmp,estimateOrder]=sort(-Iq);

% plot number of estimates in each cluster
subplot(1,2,2);
barh(s.N(estimateOrder));
set(gca,'ytick',1:Ncluster,'yticklabel', ...
	clusterlabels(estimateOrder),'ydir','reverse');
text(s.N(estimateOrder),1:Ncluster,cellstr(num2str(s.N(estimateOrder)')));
axis([0 max(s.N) 0.5 Ncluster+.5]); 
title('Number of ICA estimates in the estimate-clusters');
ylabel('Label');
xlabel('Number of estimates');


set(2,'name','Icasso: Estimate Quality');


%%%%%%%%%%%%%%% Dendrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Show dendrogram 

figure(3);
clf reset;
icassoDendrogram(sR,level);

%%%%%%%%%%%%%%% Correlation graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf reset;
if ~isempty(sR.projection.coordinates),  
  icassoGraph(sR,graphArgs{:});
else
  warning(['Projection coordinates not computed, can''t start ' ...
           'icassoGraph.']);
end

%%%%%%%%%%%%% Source plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch est 
 case 'source'
  figure(5); clf reset;
  set(5,'name','Icasso: Source Estimates (centrotypes)');
  signalplot(S(estimateOrder,:)); 
  set(gca,'yticklabel',estimateOrder); 
  ylabel('Label');
  xlabel('Sample #');
  title(['Independent components (ranked according to I_q)']);
 case 'demixing'
  figure(5); clf reset;
  set(5,'name','Icasso: Demixing Matrix Estimate');
  signalplot(W(estimateOrder,:)); 
  set(gca,'yticklabel',estimateOrder); 
  ylabel('Label');
  xlabel('Column #');
  title(['Demixing matrix rows (ranked according to' ...
	 ' I_q)']);
 case 'mixing'
  figure(5); clf reset;
  set(5,'name','Icasso: Mixing Matrix Estimate');
  signalplot(A(:,estimateOrder)'); 
  set(gca,'yticklabel',estimateOrder); 
  ylabel('Label');
  xlabel('Row #');
  view([-90 90])
  title(['Mixing matrix columns (ranked according to' ...
	 ' I_q)']);
end

%%%%%%%%%%%%% Details of Iq %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if detailedrankplot,
  figure(6); clf;
  set(6,'name','Icasso: Detaild Estimate Stability');
  icassoStability(sR,level,'plotstat');  
end

%%% subfunctions 

