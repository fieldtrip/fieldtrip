function icassoDendrogram(sR,L)
%function icassoDendrogram(sR,[L])
%
%PURPOSE 
%
%To visualize the clustering and similarities between estimates as
%a dendrogram and a corrgram.
%
%EXAMPLE OF BASIC USAGE
%
% icassoDendrogram(sR,25)
%
%INPUT 
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
%sR  (struct) Icasso result data structure; clustering has to be readily
%      computed.
%[L] (scalar) number of estimate-clusters. Default: the (reduced)
%      data dimension. 
%
%DETAILS
%
%Plots dendrogram (tree visualization) of the clustering hierarchy
%and a corrgram (the values of similarity matrix visualized as a
%gray level matrix). The entries in the corrgram are ordered
%according to the leaf ordering in the dendrogram. The clusters in
%the corrgram are indicated by red lines.  
%
%SEE ALSO
% icassoGraph 
% icassoShow
% som_dendrogram
% hcluster

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

% ver 1.2 johan 100105

%% Check that clustering info is present
if isempty(sR.cluster.partition) | isempty(sR.cluster.similarity), 
  error('Missing clustering information.');
end

% Get default cluster number
if nargin<2|isempty(L),
  L=max(icassoGet(sR,'rdim'));
  disp([sprintf('\n') 'Number of estimate-clusters not given: using (reduced)' ...
	  ' data dimension.']);
end
  
%% Check that cluster L is valid
maxCluster=size(sR.cluster.partition,1);
if L<=0 | L>maxCluster,
  error('Cluster L out of range or not specified.');
end

%% Get clustering info
i=sR.cluster.dendrogram.order;
Z=sR.cluster.dendrogram.Z;
M=icassoGet(sR,'M');

% Compute cluster borders
t=zeros(1,L+1);
t(2:end-1)=find(diff(sR.cluster.partition(L,sR.cluster.dendrogram.order)))+.5;
t(1)=0.5;
t(end)=icassoGet(sR,'M')+.5;

%% Plot dendrogram
subplot(1,2,1);
som_dendrogram(sR.cluster.dendrogram.Z);set(gca,'xtick',[]);
%y=get(gca,'ytick'); y, set(gca,'yticklabel',1-y);
lim=axis;
axis([0.5 M+.5 lim(3:4)]);
view([-90 -90]);

title(sprintf('Dendrogram (linkage strategy: %s)',sR.cluster.strategy));
ylabel('Cluster merging similarity');

%% Plot similarity matrix (entries ordered according to the dendrogram)
subplot(1,2,2);
imagesc(sR.cluster.similarity(sR.cluster.dendrogram.order, ...
			      sR.cluster.dendrogram.order));
colormap(1-gray);
set(gca,'xtick',t,'ytick',t,'xticklabel','','yticklabel','', ...
	'gridline','-','xcolor','r','ycolor','r');
grid on; 
title({'Similarities between estimates',
       '(arranged according to the dendrogram)'});
xlabel('Red lines indicate estimate-clusters');

% Make labels
r=diff(t)/2+t(1:end-1);
[labels,i]=unique(sR.cluster.partition(L, ...
				    sR.cluster.dendrogram.order));
[tmp,i]=sort(i); labels=labels(i); labels=cellstr(num2str(labels(:)));

set(text(zeros(L,1),r,labels),'horiz','right');
set(ylabel({'Estimate-cluster label',' '}),'vert','bot','color','k');
%Set colorbar
colorbar;

set(gcf,'name','Icasso: Dendrogram and Similarity Matrix');
