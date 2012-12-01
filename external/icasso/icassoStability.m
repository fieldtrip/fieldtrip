function [Iq,in_avg,ext_avg,in_min,ext_max]=icassoStability(sR,L,graphmode)
%function Iq=icassoStability(sR,L,graphmode)
%
%PURPOSE
%
%To compute and/or plot the stability (quality) indices of the ICA
%estimate-clusters. 
%
% Iq=avg(intra-cluster similarity) - avg(extra-cluster similarity)
%
%See publication Himberg et al. (2004), "Validating the independent
%components of neuroimaging time-series via clustering and
%visualization". NeuroImage, 22:3(1214-1222). Ideally, each ICA
%estimate-cluster should have Iq=1. The smaller the value is, the
%less stable (compact and isolated) the estimate-cluster is.     
%
%EXAMPLES OF BASIC USAGE
%
%   Iq=icassoStability(sR);
%
%returns the stability index for each estimate (centrotype) using
%the default number of estimate-clusters. 
%
%   Iq=icassoStability(sR,13,'plotindex');
%
%...same as previous but also plots the values (in rank order) and
%uses 13 estimate clusters instead of default. 
%
%INPUT
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% sR          (struct) Icasso result data structure; clustering
%              must be readily computed. 
% [L]         (scalar) number of clusters, default: (reduced) data
%              dimension  
% [graphmode] (string) 'none' (default) | 'plotindex' | 'plotstat'   
%
%OUTPUTS
%
% Iq  (Lx1 matrix) Iq(C) contains the index value for estimate C
% (actually estimate-cluster C)
%
%DETAILS
%
%The intra-cluster similarities for cluster C mean the mutual
%similarities between the estimates in C, and the extra-cluster
%similarities for C are the similarities between the estimates in C
%and the estimates not in C. The index is the average intra-cluster
%similarity subtracted by the average intra-cluster similarity. 
%
%In default graph mode 'none'
% The function computes the index and returns it. No graphical
% output. Iq(C) contains the index for cluster C (That is, for
% estimates sR.cluster.partition(L,:)==C.
%
%In graph mode 'plotindex'
% The function computes the quality index as in mode 'none' 
% and plots the index in the active axis (gca) as follows: 
% The X-axis shows the value of the index for clusters. The clusters
% are ranked on the Y-axis in descending order according to the
% index. The Y-axis legend shows the cluster number (the same as in
% sR.cluster.partition(L,:)). The values are indicated by black
% dots connected into each other by a dotted line. 
%
%In graph mode 'plotstat'
% As previous one but instead of the index,
% some cluster statistics are shown: 
%  Light red  indicates the average intra-cluster similarity for a cluster 
%  Light blue               average extra-cluster
%  Clear red                minimum intra-cluster 
%  Clear blue               maximum extra-cluster 
%
%SEE ALSO
% clusterquality
% icassoViz
% icassoShow

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

% ver 1.2 100501 johan

if nargin<2|isempty(L),
  L=icassoGet(sR,'rdim');
  disp([sprintf('\n') 'Number of estimate-clusters not given: using (reduced)' ...
	   ' data dimension.']);
end

if nargin<3 | isempty(graphmode),
  graphmode='none';
end

clusternumber=1:L;
partition=sR.cluster.partition(L,:);

for i=1:L,
  NofEstimates(i)=sum(partition==i);
end

% compute cluster quality index

[Iq,in_avg,ext_avg]=clusterquality('mean',sR.cluster.similarity,partition);
[tmp,in_min,ext_max]=clusterquality('minmax',sR.cluster.similarity,partition);

% Sort entries according to score

[tmp,order]=sort(-Iq); 
score=Iq(order);
in_avg=in_avg(order); 
ext_avg=ext_avg(order); 
in_min=in_min(order); 
ext_max=ext_max(order);
clusternumber=clusternumber(order);
NofEstimates=NofEstimates(order);

% make tick mark labels
for i=1:L,
  ticklabel{i}=num2str(clusternumber(i));
end

% Plot 
switch lower(graphmode)
 case 'none'
  return;
 case 'plotindex'
  cla reset;
  h_score=plot(score,1:L,'o:');
  set(h_score,'color','k','markersize',8,'markerfacecolor','k');
  set(gca,'ytick',1:L,'yticklabel',ticklabel,'ydir','reve');
  ax=axis; 
  
  if ax(1)<0,
    axis([ax(1) 1  0.5 L+.5]); 
  else
    axis([0 1  0.5 L+.5]); 
  end
  
  grid on;
  xlabel('I_q=avg(S(i)_{int})-avg(S(i)_{ext})');
  title('Stability index (I_q) for ICA estimate-clusters');
 case 'plotstat'
  cla reset;
  
  h_in_avg=barh(in_avg,0.7); hold on
  h_ext_max=barh(-ext_max,0.35); 
  h_in_min=barh(in_min,0.35); 
  h_ext_avg=barh(-ext_avg,0.7); 
  
  axis on;
  set(h_in_avg,'facecolor',[1 0.7 0.7]); hold on;
  set(h_in_min,'facecolor',[1 0 0]);
  set(h_ext_avg,'facecolor',[0.7 0.7 1]);
  set(h_ext_max,'facecolor',[0 0 1]);
  axis([-1 1 0.5 L+.5]); 
  
  set(gca,'ytick',1:L,'yticklabel',ticklabel,...
	  'ydir','reve','xtick',[-1 -.5 0 .5 1],...
	  'xgrid','on',...
	  'xticklabel',{'1' '0.5 ' '0' '0.5' '1'});
  
  title('Statistics on within- and between-cluster similarities');
  xlabel('clusters are ordered according to I_q=avg(S(i)_{int})-avg(S(i)_{ext})'); 
  legend([h_in_avg(1),h_ext_avg(1),h_in_min(1),h_ext_max(1)],...
	 {'mean S_{in}' 'mean S_{ex}' 'min S_{in}' 'max S_{ex}'},-1);
 otherwise
  error('Graphmode must be ''none'',''plotindex'' or ''plotstat''.');
end

ylabel('Label')



  

