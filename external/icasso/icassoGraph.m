function icassoGraph(sR, varargin)
%function icassoGraph(sR,['identifier1',value1,'indentifier2',value2,...])
%
%PURPOSE
%
%To visualize the estimate space as a 2D projection and the
%similarities between the estimates as a graph. Estimates are
%presented as black points that are located so that the distances
%between the points approximate the similarities between the
%estimates. Estimates belonging to the same cluster are bounded by
%a red convex hull whose background color expresses its average
%density.  
%
%EXAMPLES OF BASIC USAGE
%
%Show results for as many estimate-clusters as there are dimensions
%(possibly reduced) and use default visualization  parameters: 
%
%  icassoGraph(sR); 
%
%Show results for 9 estimate-clusters. Suppress graph lines for
%similarities under value 0.7 everywhere as well as inside
%those clusters whose intra-cluster similarity is above 0.8:   
%
%  icassoGraph(sR,'L',9,'graphlimit',0.7,'dense',0.8);
%
%INPUTS
%
%sR (struct) Icasso result structure
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)  
%
% 'L' (string) 'rdim' (default) | (integer) 
%   sets the number of estimate-clusters 'rdim' sets it equal to the
%   (reduced) data dimension 
% 'colorlimit' (vector) default [0.5 0.75 0.9] 
%   sets the thresholds for colorscale for graph lines and clusters  
% 'line' (string) 'on' (default) | 'off' 
%   whether to show the similarity graph lines or not
% 'hull' (string) 'on' (default) | 'off'
%   whether to show the "cluster hulls" or not
% 'graphlimit' (scalar) in 0...1 | (string) 'auto' (default) 
%   see additional information 2
% 'dense' (scalar) in 0...1 | (string) 'auto' (default) 
%   see additional information 3
%
%DETAILS 
%
%The function clears the active figure.
%
%Estimate-cluster labels are located at the left upper corner of
%each cluster hull. The cyan circles show the locations of the
%centrotypes (see explanation in function 'centrotype') of
%estimate-clusters on the projection. The area of a cluster is
%roughly related to its density: in general, the smaller the
%cluster area, the more compact it its and the more stable the
%corresponding IC estimate is as well. 
%
%The colorbar to the right shows the color coding of clusters and
%a floating legend for graph lines (you can drag-and-drop the
%legend if it's inconveniently placed.)
%
%The background color of each cluster hull depends on the average
%intra-cluster similarity (see additional information 1).  The color
%thresholds can be changed by input parameter 'colorlimit', e.g.,
%[0.5 0.75 0.9] sets three shades of red for 0.5...0.75 (light
%red), 0.75...0.9, and 0.9...1 (bright red). The same color coding
%is used for the graph lines that present pairwise similarities
%between the estimates. However, for sake of clarity, the
%similarity graph is, by default, suppressed inside those clusters
%whose intra-cluster similarity is above the highest threshold of
%'colorlimit' (by default 0.9). The graph is globally suppressed below
%the lowest limit (by default 0.5). This default mode can be
%changed by using parameters 'dense' and 'graphlimit' (see
%additional information 2 and 3). You can completely suppress
%completely the similarity graph completely and/or the convex
%"cluster hulls" from the figure by giving value 'off' for
%parameters 'line' and/or 'hull'. Note that the function may
%automatically change these threshold in order to keep the number
%of lines reasonable, i.e., less than 5000.       
%
%ADDITIONAL INFORMATION
%
%(1)The intra-cluster similarities for cluster C are
%the mutual similarities between the estimates in C, and the
%extra-cluster similarities for C are the similarities between the
%estimates in C and the estimates not in C.
%
%(2)Parameter 'graphlimit' sets the lowest similarity value for
%which the graph lines are drawn. The default mode (string 'auto')
%sets this limit equal to the lowest 'colorlimit (0.5 by default). A
%higher or lower limit can be given explicitly. When you set a
%'graphlimit' lower than the lowest 'colorlimit' the lines in
%between appear as thin gray lines instead of red lines.  
%
%(3)Parameter 'dense' sets the limit for intra-cluster similarity
%above which above which the similarity graph inside a cluster hull
%is suppressed. String 'auto' sets the value equal to highest
%'colorlimit' (0.9 by default).    
%
%SEE ALSO
% icassoShow
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

% ver 1.2 100105

if nargin<1|isempty(sR),
  error('At least one input argument expected');
end

if isempty(sR.projection.coordinates) | isempty(sR.cluster.partition) | ...
      isempty(sR.cluster.index) | isempty(sR.cluster.similarity), 
  error('Missing similarity/projection/cluster information.');
end

% initiate output args.
index2centrotypes=[]; clusterquality=[]; partition=[];

EDGECOLOR=[0.3 0.3 0.3];

%% Set defaults and process optional input
default={'line','on','l', 'rdim',...
	 'graphlimit','auto','colorlimit',[0.5 0.75 0.9],...
	 'dense','auto','hull','on'};

varargin=processvarargin(varargin,default);
num_of_args=length(varargin);

for i=1:2:num_of_args,
  id=varargin{i}; value=varargin{i+1};
  switch lower(id)
   case 'line'
    switch lower(value)
     case 'on'
      graph=1;
     case 'off'
      graph=0; 
     otherwise
      error('Option ''line'' must be ''on'' or ''off''.');
    end
   case 'hull'
    switch lower(value)
     case 'on'
      hull=1;
     case 'off'
      hull=0;
     otherwise
      error('Option ''hull'' must be ''on'' or ''off''.');
    end
   case 'l'
    if isnumeric(value);
      level=value;
    else
      switch lower(value)
       case 'rdim'
	level=icassoGet(sR,'rdim');
       otherwise
	error(['Unknown value for identifier ' id]);
      end
    end
   case 'dense'
    if ischar(value),
      switch lower(value)
       case 'auto'
	internallimit='auto';
       otherwise
	error('Option ''dense'' must be string ''auto'' or a scalar in 0...1');
      end
    else
      internallimit=value(1); 
      if internallimit<0 | internallimit>1,
	error('Option ''dense'' must be string ''auto'' or a scalar in 0...1');
      end
    end
   case 'graphlimit'
    if ischar(value),
      switch lower(value)
       case 'auto'
	lowlimit='auto';
       otherwise
	error('Option ''graphlimit'' must be string ''auto'' or a scalar in 0...1');
      end
    else
      lowlimit=value(1); 
      if lowlimit<0 | lowlimit>1,
	error('Option ''graphlimit'' must be string ''auto'' or a scalar in 0...1');
      end
    end
   case 'colorlimit'
    if any(value(:)==0) | any(value(:)==1),
      error('0 and 1 not allowed in ''colorlimit''.');
    end
    clustercolorlimit=value(:);
   otherwise
    error(['Doesn''t recognize option ''' id '''.' sprintf('\n')...
	   'Available: ''level'',''dense'',''graphlimit'',' ...
	   '''colorlimit'',''hull'', and ''line''.']);
  end
end

% set auto values
if strcmp(lowlimit,'auto'),
  lowlimit=min(clustercolorlimit);
end

if strcmp(internallimit,'auto'),
  internallimit=max(clustercolorlimit);
end

% Check if cluster level is valid
maxCluster=size(sR.cluster.partition,1);
if level<=0 | level>maxCluster,
  error('Cluster level out of range or not specified.');
end

%%%%% Compute some cluster statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the partition

partition=sR.cluster.partition(level,:);
Ncluster=max(partition);

% cluster statistics
c=sR.cluster.similarity;
s=clusterstat(c,partition);

%%%% Get centrotypes %

index2centrotypes=icassoIdx2Centrotype(sR,'partition',partition);

%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf reset; hold on;

p=sR.projection.coordinates; 

% define clustercolors
clustercolor=redscale(length(clustercolorlimit)+1);

% initiate graphic handles
h_graph=[];   % graph lines
graphText=[]; % label texts for graph lines

% Reduce similarities 
% If partitioning is computed, limit not only by lowlimit but also
% by denseLim: ignore lines inside cluster hulls if average internal
% similarities are over denseLim (dense clusters)

if graph,
  c=reducesim(c,lowlimit,partition,s.internal.avg,internallimit);
end

% set colors for clusters

faceColorMatrix=repmat(NaN,Ncluster,3);
for i=1:length(clustercolorlimit);
  tmp=find(s.internal.avg(:)>=clustercolorlimit(i));
  faceColorMatrix(tmp,:)=repmat(clustercolor(i+1,:),length(tmp),1);
end

% set edgecolors
edgeColorMatrix=repmat(EDGECOLOR,Ncluster,1);
  
% draw faces for clusters; they have to be in
% bottom; otherwise they shade everything else

if hull,
  h_fill=clusterhull('fill',p,partition,faceColorMatrix);
end

title(sprintf('Estimate space as a 2D %s projection',upper(sR.projection.method)));

if lowlimit<clustercolorlimit(1),
  graphlimit=[lowlimit, clustercolorlimit(:)', 1];
  linecolor=[repmat(clustercolor(2,2),1,3).^.5; clustercolor(2:end,:)];
else
  graphlimit=[clustercolorlimit(:)', 1];
  linecolor=clustercolor(2:end,:);
end

if graph,
  h_graph=similaritygraph(p,c, graphlimit, linecolor); 
  ax=axis;
  xwidth=ax(2)-ax(1);
else
  % Only vertices
  h_graph=similaritygraph(p);
  ax=axis;
  xwidth=ax(2)-ax(1);
end

% Cluster labels
txt=cellstr(num2str([1:Ncluster]'));
% Hull edges...

if hull,
  [h_edge,txtcoord]=clusterhull('edge',p,partition,edgeColorMatrix); 
  %h_clusterlabel=text(txtcoord(:,1)-xwidth/100,txtcoord(:,2),txt);
end

%...and centrotypes (cyan circles)

%% Plot centrotypes
h_centrotype=plot(p(index2centrotypes,1),p(index2centrotypes,2),'ko');
set(h_centrotype,'markersize',10,'color','c');

h_clusterlabel=text(p(index2centrotypes,1)-xwidth/100,p(index2centrotypes,2),txt);
for i=1:Ncluster,
  set(h_clusterlabel(i),'horiz','right','color',[0 0 0.7],'fontsize',17);
end

set(gca,'box','on');

%% Set colorbar that shows the cluster colors:

if hull,
  caxis([0 size(clustercolor,1)]); 
  colormap(clustercolor); 
  h_colorbar=colorbar('vert');
  set(h_colorbar,'ytick',[0:size(clustercolor,1)]',...
		 'yticklabel',[0;clustercolorlimit(:); 1],'dataaspect',[1 1 ...
		    1],'plotboxaspect',[1 1 1]);
  set(get(h_colorbar,'ylabel'),'string','Averge intra-cluster similarity (cluster compactness)');
  % Set legend
end

if graph,
  i=1;
  while i<=length(h_graph.example)
    if isnan(h_graph.example(i)),
      h_graph.example(i)=[];
      h_graph.edge(i)=[];
      h_graph.text(i)=[];
    else 
      i=i+1;
    end
  end
  legendText{1}='Single-run-estimate'; 
  legendText{2}='"Best Estimate" (centrotype)'; 
  legendText=[legendText, h_graph.text];
  legendHandle=[h_graph.vertex(1);h_centrotype(1);h_graph.example(:)];
else
  legendText{1}='Single-run-estimate'; 
  legendText{2}='"Best Estimate" (centrotype)'; 
  legendHandle=[h_graph.vertex(1);h_centrotype(1)];
end

legend(legendHandle,legendText);

% Set figure name

if graph & internallimit<1,
  xlabel(['Note that the pairwise similarity graph between estimates' ...
	  ' inside clusters' ...
	  sprintf('\n') ...
	  'is omitted if' ...
	  ' the average intra-cluster similarity is above ' sprintf('%0.2f',internallimit)]); 
else
  ;
end
if hull,
ylabel(['Convex hulls represent estimate-clusters.',...
	sprintf('\n'), 'Compact and isolated clusters suggest reliable estimates']);
end
set(gca,'yaxisloc','right');
set(gcf,'name','Icasso: Similarity Graph');
